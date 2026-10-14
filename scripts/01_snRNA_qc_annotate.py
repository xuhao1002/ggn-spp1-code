#!/usr/bin/env python
"""
第一步（上）：Peng 2026 单核转录组（GSE308103）质控、去双细胞、去批次、聚类、大类细胞注释。

用法（Mac 终端）：
    eval "$(conda shell.bash hook)"; conda activate ggn_target
    cd ~/全部数据库/GGN_target && python scripts/01_snRNA_qc_annotate.py

输出：results/01_snRNA/
    sample_meta.tsv                    每个样本的患者、阶段、文件名
    qc_per_sample.tsv                  每个样本质控前后细胞数、中位基因数/计数、双细胞比例
    qc_violin.png                      质控指标分布（按样本）
    cluster_annotation_auto.tsv        每个簇的自动大类标签、得分差、样本/患者分布（含空列 final_label 供人工复核）
    cluster_marker_scores.tsv          每个簇 × 每类细胞的标志基因平均得分
    cluster_top_markers.tsv            每个簇前 15 个标志基因（仅用于注释，不做统计推断）
    composition_counts_per_sample.tsv  每个样本各大类细胞数
    composition_frac_per_sample.tsv    每个样本各大类细胞占比
    umap_*.png, dotplot_markers.png    图
    snRNA_annotated.h5ad               含原始计数（layers['counts']）的完整对象
    run_log.txt                        参数、软件版本、各步耗时

说明：
- 数据为 10x Flex 探针法（18,082 个基因，无线粒体基因），故质控不用线粒体比例。
- 质控阈值：每个样本内按 MAD 判定离群（log 计数与 log 基因数，5 倍 MAD），并设下限
  （≥200 个基因、≥400 个计数）；scrublet 按样本分别去双细胞。
- 自动注释只给大类，最终标签需人工复核 cluster_annotation_auto.tsv 后填 final_label。
- 中间结果缓存在 results/01_snRNA/cache/，重跑时自动跳过读取步骤。
"""
import os, re, sys, glob, time, json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/全部数据库/GGN_target")
RAW = os.path.join(BASE, "raw/snRNA_GSE308103")
OUT = os.path.join(BASE, "results/01_snRNA")
CACHE = os.path.join(OUT, "cache")
os.makedirs(CACHE, exist_ok=True)

SEED = 0
N_WORKERS = min(6, os.cpu_count() or 1)
STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]

# ---------- 日志 ----------
LOG = open(os.path.join(OUT, "run_log.txt"), "a", encoding="utf-8")
T0 = time.time()
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

# ---------- 1. 样本信息 ----------
def parse_name(path):
    # GSM9237905_P4_Normal1.raw_counts.mtx.txt.gz -> GSM9237905, P4, Normal, P4_Normal1
    fn = os.path.basename(path)
    m = re.match(r"^(GSM\d+)_(P\d+)_([A-Za-z]+?)(\d*)\.raw_counts\.mtx\.txt\.gz$", fn)
    if not m:
        raise ValueError(f"无法解析文件名：{fn}")
    gsm, pt, stage, rep = m.groups()
    if stage not in STAGE_ORDER:
        raise ValueError(f"未知阶段 {stage}：{fn}")
    return dict(gsm=gsm, patient=pt, stage=stage, rep=rep or "0",
                sample=f"{pt}_{stage}{rep}", file=fn)

def read_dense(path):
    """致密文本矩阵：首行 = 细胞条形码（比数据行少 1 列），之后每行 = 基因名 + 各细胞计数。"""
    import anndata as ad
    info = parse_name(path)
    cache = os.path.join(CACHE, info["sample"] + ".h5ad")
    if os.path.exists(cache):
        return cache
    import gzip
    with gzip.open(path, "rt") as h:
        barcodes = h.readline().rstrip("\n").split("\t")
    df = pd.read_csv(path, sep="\t", skiprows=1, header=None, index_col=0, engine="c")
    if df.shape[1] != len(barcodes):
        raise ValueError(f"{info['file']}：列数 {df.shape[1]} 与条形码数 {len(barcodes)} 不符")
    X = sp.csr_matrix(df.to_numpy(dtype=np.float32).T)
    genes = pd.Index(df.index.astype(str).tolist())   # index_col=0 会把索引名设成整数 0，h5ad 不接受，故重建无名索引
    a = ad.AnnData(X=X,
                   obs=pd.DataFrame(index=pd.Index([f"{info['sample']}:{b}" for b in barcodes])),
                   var=pd.DataFrame(index=genes))
    for k in ("gsm", "patient", "stage", "sample"):
        a.obs[k] = info[k]
    tmp = cache + ".tmp"            # 先写临时文件再改名，避免中断后留下损坏的缓存
    a.write_h5ad(tmp)
    os.replace(tmp, cache)
    return cache

def main():
    import scanpy as sc
    import anndata as ad
    sc.settings.verbosity = 1
    sc.settings.n_jobs = N_WORKERS
    sc.settings.figdir = OUT
    log("===== 开始 =====")
    log(f"scanpy {sc.__version__}, anndata {ad.__version__}, numpy {np.__version__}, pandas {pd.__version__}")

    files = sorted(glob.glob(os.path.join(RAW, "GSM*.raw_counts.mtx.txt.gz")))
    meta = pd.DataFrame([parse_name(f) for f in files])
    log(f"找到 {len(files)} 个样本，{meta.patient.nunique()} 名患者；阶段：{meta.stage.value_counts().to_dict()}")
    meta.to_csv(os.path.join(OUT, "sample_meta.tsv"), sep="\t", index=False)

    # ---------- 2. 读取（并行，带缓存） ----------
    merged = os.path.join(CACHE, "merged_raw.h5ad")
    if os.path.exists(merged):
        adata = ad.read_h5ad(merged); log(f"读取缓存 merged_raw.h5ad：{adata.shape}")
    else:
        with ProcessPoolExecutor(N_WORKERS) as ex:
            paths = list(ex.map(read_dense, files))
        parts = [ad.read_h5ad(p) for p in paths]
        genes0 = parts[0].var_names
        for p, a in zip(paths, parts):
            if not a.var_names.equals(genes0):
                raise ValueError(f"{p} 的基因顺序与第一个样本不同")
        adata = ad.concat(parts, join="outer", merge="same")
        if adata.var_names.has_duplicates:
            log(f"重复基因名 {adata.var_names.duplicated().sum()} 个，已加后缀去重")
            adata.var_names_make_unique()
        adata.write_h5ad(merged)
        log(f"合并完成：{adata.n_obs} 个细胞 × {adata.n_vars} 个基因")

    for k in ("patient", "stage", "sample", "gsm"):
        adata.obs[k] = adata.obs[k].astype(str)
    adata.obs["stage"] = pd.Categorical(adata.obs["stage"], categories=STAGE_ORDER, ordered=True)

    # ---------- 3. 质控 ----------
    sc.pp.calculate_qc_metrics(adata, percent_top=[20], log1p=True, inplace=True)
    n_before = adata.obs.groupby("sample", observed=True).size()

    def mad_outlier(x, nmad=5):
        med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826
        return (x < med - nmad * mad) | (x > med + nmad * mad)

    adata.obs["qc_outlier"] = False
    for s, idx in adata.obs.groupby("sample", observed=True).indices.items():
        o = adata.obs.iloc[idx]
        bad = (mad_outlier(o["log1p_total_counts"].values)
               | mad_outlier(o["log1p_n_genes_by_counts"].values)
               | mad_outlier(o["pct_counts_in_top_20_genes"].values))
        adata.obs.iloc[idx, adata.obs.columns.get_loc("qc_outlier")] = bad
    adata.obs["qc_low"] = (adata.obs["n_genes_by_counts"] < 200) | (adata.obs["total_counts"] < 400)

    # 质控图（过滤前）
    fig, axes = plt.subplots(3, 1, figsize=(22, 10), sharex=True)
    meta["_p"] = meta["patient"].str[1:].astype(int)
    meta["_s"] = meta["stage"].map({s: i for i, s in enumerate(STAGE_ORDER)})
    order = meta.sort_values(["_p", "_s", "rep"])["sample"].tolist()
    meta.drop(columns=["_p", "_s"], inplace=True)
    for ax, col in zip(axes, ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes"]):
        data = [adata.obs.loc[adata.obs["sample"] == s, col].values for s in order]
        ax.boxplot(data, showfliers=False); ax.set_ylabel(col, fontsize=8)
    axes[-1].set_xticks(range(1, len(order) + 1)); axes[-1].set_xticklabels(order, rotation=90, fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "qc_violin.png"), dpi=130); plt.close(fig)

    keep = ~(adata.obs["qc_outlier"] | adata.obs["qc_low"])
    log(f"质控去除 {(~keep).sum()} 个细胞（低质量 {adata.obs['qc_low'].sum()}，MAD 离群 {adata.obs['qc_outlier'].sum()}，两者可重叠）")
    adata = adata[keep.values].copy()
    sc.pp.filter_genes(adata, min_cells=10)
    log(f"基因过滤后：{adata.n_vars} 个基因")

    # ---------- 4. 去双细胞（按样本） ----------
    adata.layers["counts"] = adata.X.copy()
    sc.pp.scrublet(adata, batch_key="sample", random_state=SEED)
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    dbl = adata.obs.groupby("sample", observed=True)["predicted_doublet"].mean()
    log(f"scrublet 预测双细胞 {adata.obs['predicted_doublet'].sum()} 个（各样本中位 {dbl.median():.1%}，最高 {dbl.max():.1%}）")
    n_after_qc = adata.obs.groupby("sample", observed=True).size()
    adata = adata[~adata.obs["predicted_doublet"].values].copy()

    qc = pd.DataFrame({
        "cells_raw": n_before, "cells_after_qc": n_after_qc,
        "doublet_rate": dbl, "cells_final": adata.obs.groupby("sample", observed=True).size(),
        "median_genes": adata.obs.groupby("sample", observed=True)["n_genes_by_counts"].median(),
        "median_counts": adata.obs.groupby("sample", observed=True)["total_counts"].median(),
    }).reindex(order)
    qc = meta.set_index("sample")[["patient", "stage", "gsm"]].join(qc)
    qc.to_csv(os.path.join(OUT, "qc_per_sample.tsv"), sep="\t")
    low = qc[qc["cells_final"] < 200]
    log(f"最终 {adata.n_obs} 个细胞；细胞数 < 200 的样本：{list(low.index) if len(low) else '无'}")

    # ---------- 5. 归一化、高变基因、PCA、Harmony ----------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3",
                                    layer="counts", batch_key="sample")
        hvg_flavor = "seurat_v3"
    except Exception as e:  # 未装 scikit-misc 时退回
        log(f"seurat_v3 不可用（{type(e).__name__}），改用 seurat 法")
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="sample")
        hvg_flavor = "seurat"
    log(f"高变基因 {int(adata.var['highly_variable'].sum())} 个（{hvg_flavor}，按样本）")

    ad_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(ad_hvg, max_value=10)
    sc.tl.pca(ad_hvg, n_comps=50, random_state=SEED)
    adata.obsm["X_pca"] = ad_hvg.obsm["X_pca"]
    adata.uns["pca"] = ad_hvg.uns["pca"]
    del ad_hvg

    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "sample",
                               random_state=SEED, max_iter_harmony=20)
    Z = ho.Z_corr
    if hasattr(Z, "detach"):
        Z = Z.detach().cpu().numpy()
    Z = np.asarray(Z)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    adata.obsm["X_pca_harmony"] = Z
    log(f"Harmony 完成（批次 = 样本，harmonypy {getattr(harmonypy, '__version__', '?')}）")

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    for r in (0.3, 0.6, 1.0):
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph",
                     n_iterations=2, random_state=SEED)
        log(f"leiden 分辨率 {r}：{adata.obs[f'leiden_{r}'].nunique()} 个簇")
    CL = "leiden_1.0"

    # ---------- 6. 大类注释 ----------
    MARKERS = {
        "AT2":            ["SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"],
        "AT1":            ["AGER", "RTKN2", "CLIC5", "EMP2"],
        "Club/secretory": ["SCGB1A1", "SCGB3A2", "SCGB3A1", "BPIFB1"],
        "Ciliated":       ["FOXJ1", "CAPS", "TPPP3", "RSPH1"],
        "Basal":          ["KRT5", "TP63", "KRT15"],
        "Endothelial":    ["PECAM1", "VWF", "CLDN5", "CDH5"],
        "Lymphatic EC":   ["PROX1", "CCL21", "MMRN1"],
        "Fibroblast":     ["COL1A1", "COL1A2", "DCN", "LUM"],
        "SMC/pericyte":   ["ACTA2", "MYH11", "RGS5", "PDGFRB"],
        "Mesothelial":    ["MSLN", "UPK3B", "WT1"],
        "T cell":         ["CD3E", "CD3D", "TRAC", "CD2"],
        "NK":             ["NKG7", "GNLY", "KLRD1", "KLRF1"],
        "B cell":         ["MS4A1", "CD79A", "BANK1"],
        "Plasma":         ["JCHAIN", "MZB1", "IGHG1", "XBP1"],
        "Macrophage":     ["CD68", "C1QA", "C1QB", "MARCO", "APOE"],
        "Monocyte":       ["CD14", "VCAN", "S100A8"],
        "cDC":            ["CLEC9A", "CD1C", "FCER1A", "CLEC10A"],
        "pDC":            ["LILRA4", "IL3RA", "CLEC4C"],
        "Mast":           ["CPA3", "KIT", "MS4A2", "HDC"],
        "Neutrophil":     ["CSF3R", "FCGR3B", "CXCR2"],
    }
    COMPARTMENT = {**{k: "Epithelial" for k in ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal"]},
                   **{k: "Endothelial" for k in ["Endothelial", "Lymphatic EC"]},
                   **{k: "Stromal" for k in ["Fibroblast", "SMC/pericyte", "Mesothelial"]}}
    GENERAL = ["EPCAM", "KRT8", "KRT19", "CDH1", "PTPRC", "MKI67", "TOP2A"]

    missing = {k: [g for g in v if g not in adata.var_names] for k, v in MARKERS.items()}
    missing = {k: v for k, v in missing.items() if v}
    if missing:
        log(f"标志基因不在数据中（已跳过）：{missing}")
    for k, v in MARKERS.items():
        genes = [g for g in v if g in adata.var_names]
        sc.tl.score_genes(adata, genes, score_name=f"score::{k}", random_state=SEED)
    score_cols = [f"score::{k}" for k in MARKERS]
    cs = adata.obs.groupby(CL, observed=True)[score_cols].mean()
    cs.columns = list(MARKERS)
    cs.to_csv(os.path.join(OUT, "cluster_marker_scores.tsv"), sep="\t")

    srt = np.sort(cs.values, axis=1)
    ann = pd.DataFrame({
        "auto_label": cs.idxmax(axis=1),
        "top_score": srt[:, -1].round(3),
        "second_label": cs.apply(lambda r: r.sort_values().index[-2], axis=1),
        "margin": (srt[:, -1] - srt[:, -2]).round(3),
    })
    ann["compartment"] = ann["auto_label"].map(COMPARTMENT).fillna("Immune")
    g = adata.obs.groupby(CL, observed=True)
    ann["n_cells"] = g.size()
    ann["n_samples"] = g["sample"].nunique()
    ann["n_patients"] = g["patient"].nunique()
    ann["top_patient_frac"] = g["patient"].agg(lambda x: x.value_counts(normalize=True).iloc[0]).round(3)
    st = pd.crosstab(adata.obs[CL], adata.obs["stage"], normalize="index").round(3)
    st.columns = [f"frac_{c}" for c in st.columns]
    ann = ann.join(st)
    gen = [x for x in GENERAL if x in adata.var_names]
    gx = adata[:, gen].X
    gx = gx.toarray() if sp.issparse(gx) else np.asarray(gx)
    gm = pd.DataFrame(gx, index=adata.obs_names, columns=[f"mean_{x}" for x in gen])
    ann = ann.join(gm.groupby(adata.obs[CL], observed=True).mean().round(3))
    ann["flag"] = ""
    ann.loc[ann["margin"] < 0.1, "flag"] += "得分差小;"
    ann.loc[ann["top_patient_frac"] > 0.7, "flag"] += "单一患者为主(恶性上皮候选或批次);"
    ann.loc[ann["n_cells"] < 100, "flag"] += "细胞少;"
    ann["final_label"] = ""
    ann.index.name = CL
    ann.to_csv(os.path.join(OUT, "cluster_annotation_auto.tsv"), sep="\t")
    adata.obs["auto_label"] = adata.obs[CL].map(ann["auto_label"]).astype("category")
    adata.obs["compartment"] = adata.obs[CL].map(ann["compartment"]).astype("category")
    log(f"自动注释：{adata.obs['auto_label'].value_counts().to_dict()}")

    # 簇标志基因（每簇最多抽 1000 个细胞，仅用于注释参考）
    rng = np.random.default_rng(SEED)
    idx = np.concatenate([rng.choice(v, min(len(v), 1000), replace=False)
                          for v in adata.obs.groupby(CL, observed=True).indices.values()])
    sub = adata[np.sort(idx)].copy()
    sc.tl.rank_genes_groups(sub, CL, method="wilcoxon", n_genes=15)
    top = pd.DataFrame({c: [n for n in sub.uns["rank_genes_groups"]["names"][c]]
                        for c in sub.uns["rank_genes_groups"]["names"].dtype.names}).T
    top.index.name = CL
    top.to_csv(os.path.join(OUT, "cluster_top_markers.tsv"), sep="\t")
    del sub

    # ---------- 7. 组成表 ----------
    comp = pd.crosstab(adata.obs["sample"], adata.obs["auto_label"]).reindex(order)
    comp.to_csv(os.path.join(OUT, "composition_counts_per_sample.tsv"), sep="\t")
    comp.div(comp.sum(axis=1), axis=0).round(4).to_csv(
        os.path.join(OUT, "composition_frac_per_sample.tsv"), sep="\t")

    # ---------- 8. 图 ----------
    def save_umap(color, fn, **kw):
        fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, frameon=False, **kw)
        fig.savefig(os.path.join(OUT, fn), dpi=130, bbox_inches="tight"); plt.close(fig)
    save_umap(CL, "umap_leiden1.0.png", legend_loc="on data", legend_fontsize=6)
    save_umap("auto_label", "umap_auto_label.png")
    save_umap(["stage", "patient"], "umap_stage_patient.png")
    save_umap(["log1p_total_counts", "n_genes_by_counts"], "umap_qc.png")
    genes_dot = {k: [g for g in v if g in adata.var_names][:3] for k, v in MARKERS.items()}
    dp = sc.pl.dotplot(adata, genes_dot, groupby=CL, standard_scale="var", show=False, return_fig=True)
    dp.savefig(os.path.join(OUT, "dotplot_markers.png"), dpi=130, bbox_inches="tight"); plt.close("all")

    # ---------- 9. 保存 ----------
    for c in score_cols:
        del adata.obs[c]
    adata.uns["run_params"] = json.dumps(dict(seed=SEED, hvg=hvg_flavor, harmony_key="sample",
                                              cluster_key=CL, qc="MAD5 + n_genes>=200 + counts>=400"))
    adata.write_h5ad(os.path.join(OUT, "snRNA_annotated.h5ad"), compression="gzip")
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
