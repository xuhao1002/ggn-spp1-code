#!/usr/bin/env python
"""
第二步（改进版，2026-10-08）：上皮亚群重新聚类 + infercnvpy（按患者、Normal 上皮留一参照）+ 按患者判定 CNV 恶性
+ 样本×染色体臂 CNV 谱。v1 脚本和结果保留在 scripts/02_epithelial_infercnv_v1.py、results/02_epithelial_v1/。

用法：
    cd ~/全部数据库/GGN_target && nohup /opt/anaconda3/envs/ggn_target/bin/python scripts/02_epithelial_infercnv.py > results/02_console.txt 2>&1 &

输入：results/01_snRNA/epithelial_cells.h5ad、snRNA_final.h5ad（取非上皮参照）、
      raw/annotation/gencode.v44.basic.annotation.gtf.gz、raw/annotation/cytoBand.txt.gz（hg38，着丝粒位置）
输出：results/02_epithelial/
    epi_cell_cnv.tsv                 每个上皮细胞：cnv_score、cnv_corr、cnv_malignant、cnv_high（v1 的全局阈值规则，作对照）、亚群簇
    epi_cluster_summary.tsv          上皮亚群每簇：细胞数、阶段构成、CNV 得分、cnv_malignant 比例
    epi_cnv_per_sample.tsv           每个样本：上皮细胞数、cnv_malignant 比例、cnv_high 比例（以样本为单位）
    patient_thresholds.tsv           每个患者的阈值（得分、相关）及其来源（患者自身 Normal / 全体 Normal）
    sample_arm_cnv.tsv               样本 × 染色体臂（样本平均 X_cnv 减同患者 Normal 上皮平均）
    sample_arm_cnv_vs_globalNormal.tsv  Normal 样本相对「其他患者 Normal 上皮」的臂水平偏离（查 P21_Normal）
    sample_arm_cnv_heatmap.png       行按患者分组、组内 Normal→AAH→AIS→MIA→LUAD
    cluster_cnv_heatmap.png          各上皮簇平均 X_cnv（沿基因组）
    umap_epi_*.png、cnv_score_hist.png、epi_cnv.h5ad、run_log.txt

设计（均为预先写定）：
1. 参照分两类（infercnvpy 对所有细胞用同一个「有界」规则：每个基因取两类参照均值的最小、最大值，落在区间内记 0，区间外才算偏离；
   并非逐细胞选参照，故各阶段的基准一致）：
   (a) 非上皮细胞，按 final_label 分层，每类最多 3000 个；
   (b) Normal 样本的上皮细胞（AT2、AT1、Ciliated，每类最多 3000 个），留一法：判定某患者时排除该患者自己的细胞；
       P21_Normal 疑似混入病灶细胞，不放进参照（但仍作为待判定细胞）。
2. 按患者运行 infercnv：待判定细胞 = 该患者所有上皮细胞（含其 Normal 上皮，用作该患者的零分布）。
3. 判定：该患者非 Normal 上皮中 CNV 得分前 5% 的细胞，X_cnv 取平均 = 该患者「肿瘤共识谱」；
   每个细胞与共识谱求 Pearson 相关。cnv_malignant = 得分 > 阈值 且 相关 > 阈值；
   两个阈值均为该患者 Normal 上皮的第 99 百分位；该患者 Normal 上皮 < 100 个时改用全体 Normal 上皮分布。
4. CNV 结果只作辅助指标：不进入第 1 层入选标准，也不用于切分伪批量的细胞。
"""
import os, json, time, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

BASE = os.path.expanduser("~/全部数据库/GGN_target")
IN1 = os.path.join(BASE, "results/01_snRNA")
OUT = os.path.join(BASE, "results/02_epithelial")
GTF = os.path.join(BASE, "raw/annotation/gencode.v44.basic.annotation.gtf.gz")
CYTO = os.path.join(BASE, "raw/annotation/cytoBand.txt.gz")
os.makedirs(OUT, exist_ok=True)
SEED = 0
STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]
NORMAL_REF_LABELS = ["AT2", "AT1", "Ciliated"]
EXCLUDE_REF_SAMPLES_PREFIX = "P21_Normal"
REF_MAX_PER_LABEL = 3000
TOP_FRAC = 0.05
MIN_NORMAL_CELLS = 100
CHROMS = [f"chr{i}" for i in range(1, 23)]

LOG = open(os.path.join(OUT, "run_log.txt"), "a", encoding="utf-8")
T0 = time.time()
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

def load_centromeres():
    cb = pd.read_csv(CYTO, sep="\t", header=None, names=["chrom", "start", "end", "name", "stain"])
    cen = cb[cb["stain"] == "acen"].groupby("chrom").agg(cs=("start", "min"), ce=("end", "max"))
    return ((cen["cs"] + cen["ce"]) / 2).to_dict()

def main():
    import infercnvpy as cnv
    import harmonypy
    log("===== 开始（改进版）=====")
    rng = np.random.default_rng(SEED)
    cen = load_centromeres()

    # ---------- 1. 读入 ----------
    epi = sc.read_h5ad(os.path.join(IN1, "epithelial_cells.h5ad"))
    log(f"上皮细胞：{epi.shape}")
    full = sc.read_h5ad(os.path.join(IN1, "snRNA_final.h5ad"))
    nonepi_mask = (full.obs["final_compartment"] != "Epithelial").values
    pos = np.where(nonepi_mask)[0]
    lab = full.obs["final_label"].values[pos]
    ref_idx = []
    for l in np.unique(lab):
        p = pos[lab == l]
        ref_idx.append(rng.choice(p, min(len(p), REF_MAX_PER_LABEL), replace=False))
    ref_idx = np.sort(np.concatenate(ref_idx))
    nonepi = full[ref_idx].copy()
    del full
    log(f"非上皮参照：{nonepi.n_obs} 个细胞")

    # ---------- 2. 上皮亚群重新降维聚类（仅用于分组与作图） ----------
    e2 = epi.copy()
    e2.X = e2.layers["counts"].copy()
    sc.pp.highly_variable_genes(e2, n_top_genes=3000, flavor="seurat_v3", batch_key="sample", span=0.3)
    sc.pp.normalize_total(e2, target_sum=1e4); sc.pp.log1p(e2)
    hv = e2[:, e2.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10); sc.tl.pca(hv, n_comps=50, random_state=SEED)
    epi.obsm["X_pca"] = hv.obsm["X_pca"]; del hv, e2
    ho = harmonypy.run_harmony(epi.obsm["X_pca"], epi.obs, "sample", random_state=SEED, max_iter_harmony=20)
    Z = ho.Z_corr
    Z = Z.detach().cpu().numpy() if hasattr(Z, "detach") else np.asarray(Z)
    epi.obsm["X_pca_harmony"] = Z if Z.shape[0] == epi.n_obs else Z.T
    sc.pp.neighbors(epi, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
    sc.tl.umap(epi, random_state=SEED)
    sc.tl.leiden(epi, resolution=0.5, key_added="epi_leiden", flavor="igraph", n_iterations=2,
                 directed=False, random_state=SEED)
    log(f"上皮亚群 leiden 0.5：{epi.obs['epi_leiden'].nunique()} 个簇")

    # ---------- 3. 基因坐标 ----------
    gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 8],
                      names=["chromosome", "feature", "start", "end", "attr"], low_memory=False)
    gtf = gtf[gtf["feature"] == "gene"].copy()
    gtf["gene_name"] = gtf["attr"].str.extract(r'gene_name "([^"]+)"')[0]
    gtf = gtf.dropna(subset=["gene_name"]).drop_duplicates("gene_name").set_index("gene_name")
    genes = [g for g in epi.var_names if g in gtf.index and gtf.loc[g, "chromosome"] in CHROMS and g in nonepi.var_names]
    log(f"可用基因（有坐标、chr1–22）：{len(genes)}")
    epi_g = epi[:, genes]
    non_g = nonepi[:, genes]
    var = pd.DataFrame({c: gtf.loc[genes, c].values for c in ["chromosome", "start", "end"]}, index=genes)

    X_epi = epi_g.X if sp.issparse(epi_g.X) else sp.csr_matrix(epi_g.X)
    X_non = non_g.X if sp.issparse(non_g.X) else sp.csr_matrix(non_g.X)
    obs_cols = ["sample", "patient", "stage", "final_label"]

    # ---------- 4. 按患者运行 infercnv ----------
    patients = sorted(epi.obs["patient"].unique(), key=lambda x: int(x[1:]))
    n_cells = epi.n_obs
    Xcnv_all = None
    score = pd.Series(np.nan, index=epi.obs_names)
    corr = pd.Series(np.nan, index=epi.obs_names)
    sum_by_cluster, cnt_by_cluster = {}, {}
    chr_pos = None
    cons_store, top_store, null_store = {}, {}, {}
    is_normal = (epi.obs["stage"] == "Normal").values
    sample_arr = epi.obs["sample"].astype(str).values
    patient_arr = epi.obs["patient"].astype(str).values
    label_arr = epi.obs["final_label"].astype(str).values
    nonepi_obs = nonepi.obs[obs_cols].copy()

    for p in patients:
        t0 = time.time()
        ti = np.where(patient_arr == p)[0]
        cand = np.where((~(patient_arr == p)) & is_normal & np.isin(label_arr, NORMAL_REF_LABELS)
                        & ~np.char.startswith(sample_arr.astype(str), EXCLUDE_REF_SAMPLES_PREFIX))[0]
        ri = []
        for l in NORMAL_REF_LABELS:
            c = cand[label_arr[cand] == l]
            ri.append(rng.choice(c, min(len(c), REF_MAX_PER_LABEL), replace=False))
        ri = np.concatenate(ri)
        Xc = sp.vstack([X_epi[ti], X_epi[ri], X_non]).tocsr()
        ob = pd.concat([epi.obs.iloc[ti][obs_cols], epi.obs.iloc[ri][obs_cols], nonepi_obs])
        grp = np.array(["test"] * len(ti) + ["ref_normal_epi"] * len(ri) + ["ref_nonepi"] * len(non_g))
        ad = sc.AnnData(Xc, obs=pd.DataFrame({"cnv_group": pd.Categorical(grp)}, index=[f"c{i}" for i in range(Xc.shape[0])]),
                        var=var.copy())
        cnv.tl.infercnv(ad, reference_key="cnv_group", reference_cat=["ref_normal_epi", "ref_nonepi"],
                        window_size=250, step=10, n_jobs=6)
        Xp = ad.obsm["X_cnv"]
        Xp = Xp.toarray() if sp.issparse(Xp) else np.asarray(Xp)
        Xp = Xp.astype(np.float32)[: len(ti)]
        if chr_pos is None:
            chr_pos = ad.uns["cnv"]["chr_pos"]
            Xcnv_all = np.zeros((n_cells, Xp.shape[1]), dtype=np.float32)
        Xcnv_all[ti] = Xp
        sc_p = np.abs(Xp).mean(axis=1)
        score.iloc[ti] = sc_p
        # 肿瘤共识谱：非 Normal 上皮得分前 5%
        nn = np.where(~is_normal[ti])[0]
        k = max(int(np.ceil(TOP_FRAC * len(nn))), 1)
        top = nn[np.argsort(sc_p[nn])[::-1][:k]]
        cons = Xp[top].mean(axis=0)
        cons_c = cons - cons.mean()
        Xc_c = Xp - Xp.mean(axis=1, keepdims=True)
        den = np.linalg.norm(Xc_c, axis=1) * np.linalg.norm(cons_c) + 1e-12
        corr.iloc[ti] = (Xc_c @ cons_c) / den
        nrm = np.where(is_normal[ti])[0]
        nmean = Xp[nrm].mean(axis=0)
        kn = max(int(np.ceil(TOP_FRAC * len(nrm))), 1)
        ntop = nrm[np.argsort(sc_p[nrm])[::-1][:kn]]
        cons_store[p] = cons - nmean                      # 肿瘤共识谱（相对同患者 Normal 平均）
        null_store[p] = Xp[ntop].mean(axis=0) - nmean     # 零分布：Normal 内部同法取前 5%
        top_store[p] = len(top)
        log(f"{p}：待判定 {len(ti)}（非 Normal {len(nn)}，共识谱取 {len(top)}）、参照 Normal 上皮 {len(ri)} + 非上皮 {len(non_g)}；{time.time()-t0:.0f}s")
        del Xc, ad, Xp

    epi.obs["cnv_score"] = score.values
    epi.obs["cnv_corr"] = corr.values

    # ---------- 5a. 窗口 → 染色体臂 ----------
    chrs = list(chr_pos.keys())
    starts = list(chr_pos.values()) + [Xcnv_all.shape[1]]
    gene_start = {c: np.sort(var.loc[var["chromosome"] == c, "start"].values) for c in chrs}
    win_arm = []
    for i, c in enumerate(chrs):
        nw = starts[i + 1] - starts[i]
        for j in range(nw):
            gi = int((j + 0.5) / nw * len(gene_start[c]))
            gi = min(gi, len(gene_start[c]) - 1)
            win_arm.append(f"{c[3:]}{'p' if gene_start[c][gi] < cen.get(c, np.inf) else 'q'}")
    win_arm = np.array(win_arm)
    arms = [a for a in dict.fromkeys(win_arm)]
    def arm_mean(v):
        return np.array([v[win_arm == a].mean() for a in arms])

    # ---------- 5. 阈值与判定 ----------
    ob = epi.obs
    glob_norm = ob[(ob["stage"] == "Normal")]
    g_sc, g_co = np.percentile(glob_norm["cnv_score"], 99), np.percentile(glob_norm["cnv_corr"], 99)
    rows, mal = [], pd.Series(False, index=ob.index)
    for p in patients:
        pn = ob[(ob["patient"] == p) & (ob["stage"] == "Normal")]
        if len(pn) >= MIN_NORMAL_CELLS:
            ts, tc, src = np.percentile(pn["cnv_score"], 99), np.percentile(pn["cnv_corr"], 99), "patient_normal"
        else:
            ts, tc, src = g_sc, g_co, "global_normal"
        m = (ob["patient"] == p) & (ob["cnv_score"] > ts) & (ob["cnv_corr"] > tc)
        mal[m] = True
        rows.append(dict(patient=p, n_normal_epi=len(pn), thr_score=round(ts, 5), thr_corr=round(tc, 4),
                         source=src, n_consensus_cells=top_store[p]))
    # 固定规则：共识谱的最大臂水平偏离必须超过「Normal 内部同法取前 5%」零分布的第 95 百分位（跨患者），
    # 否则该患者判为「未检出 CNV 克隆」，其全部细胞 cnv_malignant = False
    null_max = pd.Series({p: np.abs(arm_mean(null_store[p])).max() for p in patients})
    T_arm = float(np.percentile(null_max, 95))
    clone = []
    for p in patients:
        av = arm_mean(cons_store[p])
        clone.append(dict(patient=p, consensus_max_abs_arm=round(float(np.abs(av).max()), 5),
                          top_arm=arms[int(np.abs(av).argmax())], null_max_abs_arm=round(float(null_max[p]), 5),
                          arm_threshold=round(T_arm, 5), clone_detected=bool(np.abs(av).max() > T_arm),
                          n_malignant_before=int(mal[(ob["patient"] == p)].sum())))
    cl_df = pd.DataFrame(clone).set_index("patient")
    for p in cl_df.index[~cl_df["clone_detected"]]:
        mal[(ob["patient"] == p)] = False
    cl_df["n_malignant_after"] = [int(mal[(ob["patient"] == p)].sum()) for p in cl_df.index]
    cl_df.to_csv(os.path.join(OUT, "patient_cnv_clone.tsv"), sep="\t")
    log(f"臂水平阈值 {T_arm:.5f}；未检出 CNV 克隆的患者：{list(cl_df.index[~cl_df['clone_detected']]) or '无'}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "patient_thresholds.tsv"), sep="\t", index=False)
    epi.obs["cnv_malignant"] = mal.values
    # v1 规则（全体 Normal 正常上皮得分第 99 百分位）作对照
    normal_epi = ob[(ob["stage"] == "Normal") & ob["final_label"].isin(["AT2", "AT1", "Ciliated", "Basal"])]
    v1_thr = np.percentile(normal_epi["cnv_score"], 99)
    epi.obs["cnv_high"] = epi.obs["cnv_score"] > v1_thr
    ct = pd.crosstab(epi.obs["cnv_high"], epi.obs["cnv_malignant"])
    log(f"cnv_malignant 比例 {epi.obs['cnv_malignant'].mean():.3f}；cnv_high（全局阈值 {v1_thr:.5f}）比例 {epi.obs['cnv_high'].mean():.3f}")
    log("cnv_high × cnv_malignant 交叉表：" + json.dumps({f"{a}|{b}": int(v) for (a, b), v in ct.stack().items()}))
    both = ct.loc[True, True] if (True in ct.index and True in ct.columns) else 0
    log(f"一致性：两者都为真 {int(both)}；Jaccard={both / max(((epi.obs['cnv_high'] | epi.obs['cnv_malignant']).sum()), 1):.3f}")

    for col, fn in [("stage", "check_malignant_by_stage.tsv"), ("final_label", "check_malignant_by_label.tsv"),
                    ("epi_leiden", "check_malignant_by_cluster.tsv")]:
        epi.obs.groupby(col, observed=True).agg(n=("cnv_malignant", "size"), cnv_malignant=("cnv_malignant", "mean"),
                                                cnv_high=("cnv_high", "mean")).round(4).to_csv(os.path.join(OUT, fn), sep="\t")
    log("cnv_malignant 按阶段：" + json.dumps(epi.obs.groupby("stage", observed=True)["cnv_malignant"].mean().round(4).to_dict()))
    epi.obs[["sample", "patient", "stage", "final_label", "epi_leiden", "cnv_score", "cnv_corr",
             "cnv_malignant", "cnv_high"]].to_csv(os.path.join(OUT, "epi_cell_cnv.tsv"), sep="\t")
    g = epi.obs.groupby("epi_leiden", observed=True)
    summ = pd.DataFrame({
        "n_cells": g.size(), "n_patients": g["patient"].nunique(),
        "top_patient_frac": g["patient"].agg(lambda x: x.value_counts(normalize=True).iloc[0]).round(3),
        "main_final_label": g["final_label"].agg(lambda x: x.value_counts().index[0]),
        "median_cnv_score": g["cnv_score"].median().round(5),
        "cnv_malignant_frac": g["cnv_malignant"].mean().round(3), "cnv_high_frac": g["cnv_high"].mean().round(3),
    }).join(pd.crosstab(epi.obs["epi_leiden"], epi.obs["stage"], normalize="index").round(3).add_prefix("frac_"))
    summ.to_csv(os.path.join(OUT, "epi_cluster_summary.tsv"), sep="\t")
    # 各阶段可比性检查：待判定细胞中 X_cnv 恰为 0 的窗口比例（落在参照区间内或被噪声阈值置零）
    zf = pd.Series((Xcnv_all == 0).mean(axis=1), index=epi.obs_names)
    st = pd.DataFrame({"stage": epi.obs["stage"].values, "cnv_score": epi.obs["cnv_score"].values, "zero_window_frac": zf.values})
    stt = st.groupby("stage", observed=True).agg(n_cells=("cnv_score", "size"), median_cnv_score=("cnv_score", "median"),
                                                  mean_zero_window_frac=("zero_window_frac", "mean")).reindex(STAGE_ORDER).round(4)
    stt.to_csv(os.path.join(OUT, "stage_reference_comparability.tsv"), sep="\t")
    log("各阶段 CNV 得分中位数与零窗口比例：" + json.dumps(stt.to_dict("index"), ensure_ascii=False))
    ps = epi.obs.groupby("sample", observed=True).agg(
        patient=("patient", "first"), stage=("stage", "first"), n_epi=("cnv_score", "size"),
        cnv_malignant_frac=("cnv_malignant", "mean"), cnv_high_frac=("cnv_high", "mean")).round(4)
    ps = ps.assign(_o=ps["stage"].map({k: i for i, k in enumerate(STAGE_ORDER)})).sort_values(
        ["patient", "_o"], key=lambda s: s.map(lambda x: int(x[1:])) if s.name == "patient" else s).drop(columns="_o")
    ps.to_csv(os.path.join(OUT, "epi_cnv_per_sample.tsv"), sep="\t")

    # ---------- 6. 染色体臂水平的样本 CNV 谱 ----------
    patient_normal_mean = {}
    for p in patients:
        idx = np.where((patient_arr == p) & is_normal)[0]
        patient_normal_mean[p] = Xcnv_all[idx].mean(axis=0)
    prof = {}
    for s in sorted(epi.obs["sample"].astype(str).unique()):
        idx = np.where(sample_arr == s)[0]
        p = patient_arr[idx[0]]
        prof[s] = arm_mean(Xcnv_all[idx].mean(axis=0) - patient_normal_mean[p])
    arm_df = pd.DataFrame(prof, index=arms).T
    arm_df.insert(0, "n_epi", [int((sample_arr == s).sum()) for s in arm_df.index])
    arm_df.insert(0, "stage", [ob.loc[sample_arr == s, "stage"].iloc[0] for s in arm_df.index])
    arm_df.insert(0, "patient", [s.split("_")[0] for s in arm_df.index])
    arm_df.round(5).to_csv(os.path.join(OUT, "sample_arm_cnv.tsv"), sep="\t")
    # Normal 样本相对「其他患者 Normal 上皮平均」
    pn = {p: patient_normal_mean[p] for p in patients}
    rows = {}
    for s in [x for x in arm_df.index if arm_df.loc[x, "stage"] == "Normal"]:
        p = s.split("_")[0]
        others = np.mean([pn[q] for q in patients if q != p and q != "P21"], axis=0)
        idx = np.where(sample_arr == s)[0]
        rows[s] = arm_mean(Xcnv_all[idx].mean(axis=0) - others)
    ng = pd.DataFrame(rows, index=arms).T
    ng["max_abs_arm"] = ng.abs().max(axis=1)
    ng.round(5).to_csv(os.path.join(OUT, "sample_arm_cnv_vs_globalNormal.tsv"), sep="\t")
    log("Normal 样本相对其他患者 Normal 的最大臂偏离（前 5）：" +
        json.dumps(ng["max_abs_arm"].sort_values(ascending=False).head(5).round(4).to_dict()))

    # ---------- 7. 图 ----------
    def save(fig, fn):
        fig.savefig(os.path.join(OUT, fn), dpi=130, bbox_inches="tight"); plt.close(fig)
    for col, fn in [("epi_leiden", "umap_epi_leiden.png"), ("final_label", "umap_epi_label.png"),
                    ("cnv_score", "umap_epi_cnv_score.png"), ("cnv_malignant", "umap_epi_cnv_malignant.png"),
                    ("stage", "umap_epi_stage.png"), ("patient", "umap_epi_patient.png")]:
        save(sc.pl.umap(epi, color=col, show=False, return_fig=True, frameon=False,
                        legend_loc="on data" if col == "epi_leiden" else "right margin"), fn)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    nm = (ob["stage"] == "Normal").values
    ax.hist(ob.loc[nm, "cnv_score"], bins=100, alpha=.5, density=True, label="Normal epithelium")
    ax.hist(ob.loc[~nm, "cnv_score"], bins=100, alpha=.4, density=True, label="non-Normal epithelium")
    ax.axvline(v1_thr, color="k", ls="--"); ax.legend(); ax.set_xlabel("CNV score")
    save(fig, "cnv_score_hist.png")
    # 样本臂热图
    pm = arm_df[arm_df["stage"] != "Normal"].copy()
    pm["_o"] = pm["stage"].map({k: i for i, k in enumerate(STAGE_ORDER)})
    pm["_p"] = pm["patient"].map(lambda x: int(x[1:]))
    pm = pm.sort_values(["_p", "_o"])
    M = pm[arms].values
    vmax = np.nanpercentile(np.abs(M), 99)
    fig, ax = plt.subplots(figsize=(14, max(6, 0.17 * len(pm))))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(arms))); ax.set_xticklabels(arms, rotation=90, fontsize=7)
    ax.set_yticks(range(len(pm))); ax.set_yticklabels([f"{i} ({s})" for i, s in zip(pm.index, pm["stage"])], fontsize=6)
    prev = None
    for k, pt in enumerate(pm["patient"]):
        if prev is not None and pt != prev:
            ax.axhline(k - 0.5, color="k", lw=0.6)
        prev = pt
    plt.colorbar(im, ax=ax, shrink=0.5, label="CNV 相对同患者 Normal 上皮")
    save(fig, "sample_arm_cnv_heatmap.png")
    # 簇热图（沿基因组）
    cl = epi.obs["epi_leiden"].astype(str).values
    ucl = sorted(set(cl), key=int)
    CM = np.vstack([Xcnv_all[cl == c].mean(axis=0) for c in ucl])
    fig, ax = plt.subplots(figsize=(14, 4))
    v = np.nanpercentile(np.abs(CM), 99)
    ax.imshow(CM, aspect="auto", cmap="RdBu_r", vmin=-v, vmax=v)
    ax.set_yticks(range(len(ucl))); ax.set_yticklabels(ucl)
    ax.set_xticks([s + (e - s) / 2 for s, e in zip(starts[:-1], starts[1:])])
    ax.set_xticklabels([c[3:] for c in chrs], fontsize=7)
    for s in starts: ax.axvline(s, color="k", lw=0.3)
    save(fig, "cluster_cnv_heatmap.png")

    # ---------- 8. 保存 ----------
    epi.obsm["X_cnv"] = Xcnv_all
    epi.uns["cnv_params"] = json.dumps(dict(seed=SEED, window=250, step=10, ref_max_per_label=REF_MAX_PER_LABEL,
                                            top_frac=TOP_FRAC, v1_thr=float(v1_thr)))
    epi.write_h5ad(os.path.join(OUT, "epi_cnv.h5ad"), compression="lzf")
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
