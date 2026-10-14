#!/usr/bin/env python
"""
第二步：应用确认后的注释 → 剔除双细胞/血小板簇 → 上皮亚群重新降维聚类 → infercnvpy 推断拷贝数 → 初步区分恶性上皮。

用法：
    cd ~/全部数据库/GGN_target && /opt/anaconda3/envs/ggn_target/bin/python scripts/02_epithelial_infercnv.py

输入：results/01_snRNA/snRNA_annotated.h5ad、cluster_annotation_auto.tsv（final_label 列）、
      raw/annotation/gencode.v44.basic.annotation.gtf.gz
输出：results/02_epithelial/
    cell_labels.tsv                   每个细胞的 sample/patient/stage/final_label（剔除后）
    composition_final_per_sample.tsv  每个样本各 final_label 细胞数、占比
    epi_cluster_summary.tsv           上皮亚群每个簇：细胞数、阶段构成、CNV 得分、CNV 高细胞比例
    epi_cnv_per_sample.tsv            每个样本上皮细胞数、CNV 高细胞比例（以样本为单位）
    epi_cell_cnv.tsv                  每个上皮细胞的 CNV 得分与亚群簇
    umap_epi_*.png, cnv_heatmap.png, cnv_score_hist.png
    epi_cnv.h5ad                      上皮细胞对象（含 X_cnv、umap、cnv_score）
    run_log.txt

恶性判定（初步，预先写死）：以 Normal 阶段样本中上皮细胞 CNV 得分的第 99 百分位为阈值，
超过者记为「CNV 高」。最终是否恶性需结合簇层面比例和人工核对，不据单细胞阈值直接下结论。
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

BASE = os.path.expanduser("~/全部数据库/GGN_target")
IN1 = os.path.join(BASE, "results/01_snRNA")
OUT = os.path.join(BASE, "results/02_epithelial")
GTF = os.path.join(BASE, "raw/annotation/gencode.v44.basic.annotation.gtf.gz")
os.makedirs(OUT, exist_ok=True)
SEED = 0
CL = "leiden_1.0"
STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]
DROP = ["Doublet/ambiguous", "Platelet/Megakaryocyte"]
EPI = ["AT2", "AT1", "Ciliated", "Basal", "Proliferating epithelial",
       "AT2/Club-like epithelial", "AT2-like (patient-specific)"]
NORMAL_EPI = ["AT2", "AT1", "Ciliated", "Basal"]
REF_MAX_PER_LABEL = 3000

LOG = open(os.path.join(OUT, "run_log.txt"), "a", encoding="utf-8")
T0 = time.time()
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

def main():
    import infercnvpy as cnv
    log("===== 开始 =====")
    log(f"scanpy {sc.__version__}, infercnvpy {getattr(cnv, '__version__', '?')}")
    adata = sc.read_h5ad(os.path.join(IN1, "snRNA_annotated.h5ad"))
    log(f"读取：{adata.shape}")

    # ---------- 1. 应用 final_label ----------
    ann = pd.read_csv(os.path.join(IN1, "cluster_annotation_auto.tsv"), sep="\t", index_col=0)
    assert ann["final_label"].notna().all() and (ann["final_label"] != "").all(), "final_label 有空值"
    ann.index = ann.index.astype(str)
    adata.obs["final_label"] = adata.obs[CL].astype(str).map(ann["final_label"]).astype("category")
    n0 = adata.n_obs
    adata = adata[~adata.obs["final_label"].isin(DROP)].copy()
    log(f"剔除 {DROP}：{n0} → {adata.n_obs} 个细胞")
    adata.obs["final_label"] = adata.obs["final_label"].astype(str)
    lab = adata.obs[["sample", "patient", "stage", CL, "final_label"]]
    lab.to_csv(os.path.join(OUT, "cell_labels.tsv"), sep="\t")
    comp = pd.crosstab(lab["sample"], lab["final_label"])
    comp.to_csv(os.path.join(OUT, "composition_final_per_sample.tsv"), sep="\t")
    comp.div(comp.sum(axis=1), axis=0).round(4).to_csv(
        os.path.join(OUT, "composition_final_frac_per_sample.tsv"), sep="\t")
    log(f"final_label 计数：{lab['final_label'].value_counts().to_dict()}")

    # ---------- 2. 参照细胞（非上皮，分层抽样） ----------
    rng = np.random.default_rng(SEED)
    is_epi = adata.obs["final_label"].isin(EPI).values
    ref_idx = []
    for l, idx in adata.obs[~is_epi].groupby("final_label").indices.items():
        pos = np.where(~is_epi)[0][idx]
        ref_idx.append(rng.choice(pos, min(len(pos), REF_MAX_PER_LABEL), replace=False))
    ref_idx = np.concatenate(ref_idx)
    log(f"上皮 {is_epi.sum()} 个；参照（非上皮抽样）{len(ref_idx)} 个")

    # ---------- 3. 上皮亚群重新降维聚类 ----------
    epi = adata[is_epi].copy()
    epi.X = epi.layers["counts"].copy()
    sc.pp.highly_variable_genes(epi, n_top_genes=3000, flavor="seurat_v3", batch_key="sample", span=0.3)
    sc.pp.normalize_total(epi, target_sum=1e4)
    sc.pp.log1p(epi)
    hv = epi[:, epi.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=50, random_state=SEED)
    epi.obsm["X_pca"] = hv.obsm["X_pca"]
    del hv
    import harmonypy
    ho = harmonypy.run_harmony(epi.obsm["X_pca"], epi.obs, "sample", random_state=SEED, max_iter_harmony=20)
    Z = ho.Z_corr
    if hasattr(Z, "detach"):
        Z = Z.detach().cpu().numpy()
    Z = np.asarray(Z)
    if Z.shape[0] != epi.n_obs:
        Z = Z.T
    epi.obsm["X_pca_harmony"] = Z
    sc.pp.neighbors(epi, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
    sc.tl.umap(epi, random_state=SEED)
    sc.tl.leiden(epi, resolution=0.5, key_added="epi_leiden", flavor="igraph",
                 n_iterations=2, directed=False, random_state=SEED)
    log(f"上皮亚群 leiden 0.5：{epi.obs['epi_leiden'].nunique()} 个簇")

    # ---------- 4. infercnvpy ----------
    ref = adata[ref_idx].copy()
    ref.X = ref.layers["counts"].copy()
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    epi_c = epi.copy()
    epi_c.obs["cnv_group"] = "epithelial"
    ref.obs["cnv_group"] = "reference"
    both = sc.concat([epi_c, ref], join="inner", index_unique=None)
    both.obs["cnv_group"] = both.obs["cnv_group"].astype("category")
    del ref, epi_c
    gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 8],
                      names=["chromosome", "feature", "start", "end", "attr"], low_memory=False)
    gtf = gtf[gtf["feature"] == "gene"].copy()
    gtf["gene_name"] = gtf["attr"].str.extract(r'gene_name "([^"]+)"')[0]
    gtf = gtf.dropna(subset=["gene_name"]).drop_duplicates("gene_name").set_index("gene_name")
    for c in ["chromosome", "start", "end"]:
        both.var[c] = gtf[c].reindex(both.var_names).values
    ok = both.var["chromosome"].isin([f"chr{i}" for i in range(1, 23)])
    log(f"有染色体坐标且在 chr1–22 的基因：{int(ok.sum())} / {both.n_vars}")
    both = both[:, ok.values].copy()
    cnv.tl.infercnv(both, reference_key="cnv_group", reference_cat=["reference"],
                    window_size=250, step=10, n_jobs=6)
    log("infercnv 完成")

    # ---------- 5. CNV 得分与初步判定 ----------
    X = both.obsm["X_cnv"]
    score = np.asarray(abs(X).mean(axis=1)).ravel() if sp.issparse(X) else np.abs(X).mean(axis=1)
    both.obs["cnv_score"] = score
    is_e = (both.obs["cnv_group"] == "epithelial").values
    ref_score = score[~is_e]
    epi_score = pd.Series(score[is_e], index=both.obs_names[is_e])
    epi.obs["cnv_score"] = epi_score.reindex(epi.obs_names).values
    normal_mask = (epi.obs["stage"] == "Normal") & epi.obs["final_label"].isin(NORMAL_EPI)
    thr = float(np.percentile(epi.obs.loc[normal_mask, "cnv_score"], 99))
    epi.obs["cnv_high"] = epi.obs["cnv_score"] > thr
    log(f"CNV 阈值（Normal 正常上皮第 99 百分位）= {thr:.5f}；参照细胞中位 {np.median(ref_score):.5f}；"
        f"全部上皮中 CNV 高比例 {epi.obs['cnv_high'].mean():.3f}")
    epi.obsm["X_cnv"] = X[is_e] if not sp.issparse(X) else X[np.where(is_e)[0]]
    epi.uns["cnv"] = both.uns["cnv"]

    epi.obs[["sample", "patient", "stage", CL, "final_label", "epi_leiden", "cnv_score", "cnv_high"]].to_csv(
        os.path.join(OUT, "epi_cell_cnv.tsv"), sep="\t")

    g = epi.obs.groupby("epi_leiden", observed=True)
    summ = pd.DataFrame({
        "n_cells": g.size(), "n_patients": g["patient"].nunique(),
        "top_patient_frac": g["patient"].agg(lambda x: x.value_counts(normalize=True).iloc[0]).round(3),
        "main_final_label": g["final_label"].agg(lambda x: x.value_counts().index[0]),
        "median_cnv_score": g["cnv_score"].median().round(5),
        "cnv_high_frac": g["cnv_high"].mean().round(3),
    }).join(pd.crosstab(epi.obs["epi_leiden"], epi.obs["stage"], normalize="index").round(3).add_prefix("frac_"))
    summ.to_csv(os.path.join(OUT, "epi_cluster_summary.tsv"), sep="\t")

    ps = epi.obs.groupby("sample", observed=True).agg(
        patient=("patient", "first"), stage=("stage", "first"),
        n_epi=("cnv_score", "size"), cnv_high_frac=("cnv_high", "mean")).round(4)
    ps = ps.assign(_o=ps["stage"].map({k: i for i, k in enumerate(STAGE_ORDER)})).sort_values(["_o", "patient"]).drop(columns="_o")
    ps.to_csv(os.path.join(OUT, "epi_cnv_per_sample.tsv"), sep="\t")

    # ---------- 6. 图 ----------
    def save(fig, fn):
        fig.savefig(os.path.join(OUT, fn), dpi=130, bbox_inches="tight"); plt.close(fig)
    for col, fn in [("epi_leiden", "umap_epi_leiden.png"), ("final_label", "umap_epi_label.png"),
                    ("cnv_score", "umap_epi_cnv_score.png"), ("cnv_high", "umap_epi_cnv_high.png"),
                    ("stage", "umap_epi_stage.png"), ("patient", "umap_epi_patient.png")]:
        save(sc.pl.umap(epi, color=col, show=False, return_fig=True, frameon=False,
                        legend_loc="on data" if col == "epi_leiden" else "right margin"), fn)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(ref_score, bins=100, alpha=.5, density=True, label="reference")
    ax.hist(epi.obs.loc[normal_mask, "cnv_score"], bins=100, alpha=.5, density=True, label="Normal epithelium")
    ax.hist(epi.obs.loc[~normal_mask, "cnv_score"], bins=100, alpha=.4, density=True, label="other epithelium")
    ax.axvline(thr, color="k", ls="--"); ax.legend(); ax.set_xlabel("CNV score")
    save(fig, "cnv_score_hist.png")
    try:
        sub = both[np.sort(np.concatenate([
            rng.choice(np.where(is_e)[0], min(20000, is_e.sum()), replace=False),
            np.where(~is_e)[0][:5000]]))].copy()
        sub.obs["epi_leiden"] = epi.obs["epi_leiden"].reindex(sub.obs_names).astype(str).where(
            sub.obs["cnv_group"] == "epithelial", "reference")
        cnv.pl.chromosome_heatmap(sub, groupby="epi_leiden", show=False)
        plt.savefig(os.path.join(OUT, "cnv_heatmap.png"), dpi=110, bbox_inches="tight"); plt.close("all")
    except Exception as e:
        log(f"热图失败（不影响结果）：{e!r}")

    # ---------- 7. 保存 ----------
    epi.uns["run_params"] = json.dumps(dict(seed=SEED, cnv_thr=thr, window=250, step=10,
                                            ref_max_per_label=REF_MAX_PER_LABEL))
    epi.write_h5ad(os.path.join(OUT, "epi_cnv.h5ad"), compression="gzip")
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
