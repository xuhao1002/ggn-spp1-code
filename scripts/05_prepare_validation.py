#!/usr/bin/env python
"""
外部验证集预处理：读入 → 质控 → scrublet → 归一化/HVG/PCA/Harmony(按样本)/leiden → 标志基因自动注释。
注释沿用 01 脚本的 MARKERS（只给大类）。输出供 06 来源比值核对与外部验证使用。

用法：
    nohup /opt/anaconda3/envs/ggn_target/bin/python scripts/05_prepare_validation.py GSE164789 > results/05_GSE164789_console.txt 2>&1 &
    nohup /opt/anaconda3/envs/ggn_target/bin/python scripts/05_prepare_validation.py GSE189357 > results/05_GSE189357_console.txt 2>&1 &

输出 results/05_validation/<数据集>/：sample_meta.tsv、<数据集>_annotated.h5ad（含 layers['counts']、X=log1p）、
    cluster_annotation_auto.tsv、cluster_marker_scores.tsv、cluster_top_markers.tsv、umap_*.png、run_log.txt
质控阈值：n_genes>=200、counts>=500、线粒体<20%（10x 3' 数据有线粒体基因），scrublet 按样本去双细胞。
"""
import os, sys, re, json, time, gzip, warnings
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
warnings.filterwarnings("ignore")

DS = sys.argv[1]
BASE = os.path.expanduser("~/全部数据库/GGN_target")
OUT = os.path.join(BASE, "results/05_validation", DS)
os.makedirs(OUT, exist_ok=True)
SEED = 0
LOG = open(os.path.join(OUT, "run_log.txt"), "a", encoding="utf-8")
T0 = time.time()
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

# ---------- 样本信息 ----------
def sample_table():
    if DS == "GSE164789":
        raw = os.path.join(BASE, "raw/scRNA_GSE164789_precursor")
        meta = pd.read_csv(os.path.join(BASE, "results/05_validation/GSE164789_gsm_meta.tsv"), sep="\t")
        meta = meta[~meta["title"].str.contains("TCR")].copy()
        rows = []
        for _, r in meta.iterrows():
            pre = os.path.join(raw, f"{r['gsm']}_{r['title']}")
            rows.append(dict(sample=r["title"], gsm=r["gsm"], patient=r["title"].split("-")[0],
                             tissue="T" if r["tissue"] == "tumor" else "A", histology=r["histology"],
                             sex=r["sex"], mtx=pre + ".matrix.mtx.gz", bc=pre + ".barcodes.tsv.gz",
                             genes=pre + ".genes.tsv.gz"))
        return pd.DataFrame(rows)
    raw = os.path.join(BASE, "raw/scRNA_GSE189357")
    hist = {"TD1": "IAC", "TD2": "IAC", "TD3": "MIA", "TD4": "MIA", "TD5": "AIS", "TD6": "MIA",
            "TD7": "AIS", "TD8": "AIS", "TD9": "IAC"}
    gsm0 = 5699777
    rows = []
    for i in range(1, 10):
        t = f"TD{i}"; gsm = f"GSM{gsm0 + i - 1}"
        pre = os.path.join(raw, f"{gsm}_{t}")
        rows.append(dict(sample=t, gsm=gsm, patient=t, tissue="T", histology=hist[t], sex="",
                         mtx=pre + "_matrix.mtx.gz", bc=pre + "_barcodes.tsv.gz", genes=pre + "_features.tsv.gz"))
    return pd.DataFrame(rows)

def read_sample(args):
    name, mtx, bc, genes = args
    g = pd.read_csv(genes, sep="\t", header=None, compression="gzip")
    sym = g[1].astype(str).values
    b = pd.read_csv(bc, sep="\t", header=None, compression="gzip")[0].astype(str).values
    m = pd.read_csv(mtx, sep=r"\s+", header=None, comment="%", compression="gzip", dtype=np.int64)
    n_g, n_c, nnz = m.iloc[0].values
    m = m.iloc[1:]
    X = sp.csr_matrix((m[2].values.astype(np.float32), (m[1].values - 1, m[0].values - 1)), shape=(n_c, n_g))
    a = ad.AnnData(X, obs=pd.DataFrame(index=[f"{name}:{x}" for x in b]), var=pd.DataFrame(index=sym))
    a.var_names_make_unique()
    return a

MARKERS = {
    "AT2": ["SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"],
    "AT1": ["AGER", "RTKN2", "CLIC5", "EMP2"],
    "Club/secretory": ["SCGB1A1", "SCGB3A2", "SCGB3A1", "BPIFB1"],
    "Ciliated": ["FOXJ1", "CAPS", "TPPP3", "RSPH1"],
    "Basal": ["KRT5", "TP63", "KRT15"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Lymphatic EC": ["PROX1", "CCL21", "MMRN1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "SMC/pericyte": ["ACTA2", "MYH11", "RGS5", "PDGFRB"],
    "Mesothelial": ["MSLN", "UPK3B", "WT1"],
    "T cell": ["CD3E", "CD3D", "TRAC", "CD2"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "B cell": ["MS4A1", "CD79A", "BANK1"],
    "Plasma": ["JCHAIN", "MZB1", "IGHG1"],
    "Macrophage": ["CD68", "C1QA", "C1QB", "MARCO", "APOE"],
    "Monocyte": ["CD14", "VCAN", "S100A8"],
    "cDC": ["CLEC9A", "CD1C", "FCER1A", "CLEC10A"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C"],
    "Mast": ["CPA3", "KIT", "MS4A2", "HDC"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CXCR2"],
    "Platelet": ["PPBP", "PF4", "ITGA2B", "TUBB1"],
}
COMPARTMENT = {**{k: "Epithelial" for k in ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal"]},
               **{k: "Endothelial" for k in ["Endothelial", "Lymphatic EC"]},
               **{k: "Stromal" for k in ["Fibroblast", "SMC/pericyte", "Mesothelial"]}}

def main():
    log(f"===== 开始 {DS} =====")
    st = sample_table()
    st.drop(columns=["mtx", "bc", "genes"]).to_csv(os.path.join(OUT, "sample_meta.tsv"), sep="\t", index=False)
    log(f"样本 {len(st)} 个；{st.groupby(['tissue','histology']).size().to_dict()}")
    with ProcessPoolExecutor(max_workers=6) as ex:
        parts = list(ex.map(read_sample, [(r["sample"], r["mtx"], r["bc"], r["genes"]) for _, r in st.iterrows()]))
    adata = ad.concat(parts, join="outer", fill_value=0)
    adata.X = sp.csr_matrix(adata.X)
    del parts
    adata.obs["sample"] = adata.obs.index.str.split(":").str[0].values
    for c in ["patient", "tissue", "histology", "sex"]:
        adata.obs[c] = adata.obs["sample"].map(st.set_index("sample")[c]).values
    log(f"合并：{adata.shape}")

    # ---------- 质控 ----------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    keep = ((adata.obs["n_genes_by_counts"] >= 200) & (adata.obs["total_counts"] >= 500)
            & (adata.obs["pct_counts_mt"] < 20))
    log(f"质控去除 {int((~keep).sum())} 个细胞（{(~keep).mean()*100:.1f}%）")
    adata = adata[keep.values].copy()
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.scrublet(adata, batch_key="sample", random_state=SEED)
    nd = int(adata.obs["predicted_doublet"].sum())
    log(f"scrublet 预测双细胞 {nd} 个（{nd/adata.n_obs*100:.2f}%）")
    adata = adata[~adata.obs["predicted_doublet"].values].copy()
    log(f"质控后：{adata.shape}")
    adata.layers["counts"] = adata.X.copy()

    # ---------- 降维聚类 ----------
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", batch_key="sample", span=0.3)
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    hv = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10); sc.tl.pca(hv, n_comps=50, random_state=SEED)
    adata.obsm["X_pca"] = hv.obsm["X_pca"]; del hv
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "sample", random_state=SEED, max_iter_harmony=20)
    Z = ho.Z_corr
    Z = Z.detach().cpu().numpy() if hasattr(Z, "detach") else np.asarray(Z)
    adata.obsm["X_pca_harmony"] = Z if Z.shape[0] == adata.n_obs else Z.T
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden_1.0", flavor="igraph", n_iterations=2,
                 directed=False, random_state=SEED)
    log(f"leiden 1.0：{adata.obs['leiden_1.0'].nunique()} 个簇")

    # ---------- 自动注释 ----------
    missing = {k: [g for g in v if g not in adata.var_names] for k, v in MARKERS.items()}
    log("不在数据中的标志基因：" + json.dumps({k: v for k, v in missing.items() if v}, ensure_ascii=False))
    for k, v in MARKERS.items():
        sc.tl.score_genes(adata, [g for g in v if g in adata.var_names], score_name=f"score::{k}", random_state=SEED)
    sc_cols = [f"score::{k}" for k in MARKERS]
    CL = "leiden_1.0"
    cs = adata.obs.groupby(CL, observed=True)[sc_cols].mean(); cs.columns = list(MARKERS)
    cs.to_csv(os.path.join(OUT, "cluster_marker_scores.tsv"), sep="\t")
    srt = np.sort(cs.values, axis=1)
    ann = pd.DataFrame({"auto_label": cs.idxmax(axis=1), "top_score": srt[:, -1].round(3),
                        "second_label": cs.apply(lambda r: r.sort_values().index[-2], axis=1),
                        "margin": (srt[:, -1] - srt[:, -2]).round(3)})
    ann["compartment"] = ann["auto_label"].map(COMPARTMENT).fillna("Immune")
    g = adata.obs.groupby(CL, observed=True)
    ann["n_cells"] = g.size(); ann["n_samples"] = g["sample"].nunique(); ann["n_patients"] = g["patient"].nunique()
    ann["top_patient_frac"] = g["patient"].agg(lambda x: x.value_counts(normalize=True).iloc[0]).round(3)
    gen = [x for x in ["EPCAM", "KRT8", "KRT19", "CDH1", "PTPRC", "MKI67", "TOP2A"] if x in adata.var_names]
    gx = adata[:, gen].X; gx = gx.toarray() if sp.issparse(gx) else np.asarray(gx)
    gm = pd.DataFrame(gx, index=adata.obs_names, columns=[f"mean_{x}" for x in gen])
    ann = ann.join(gm.groupby(adata.obs[CL].values).mean().round(3))
    ann.index.name = CL
    ann["final_label"] = ann["auto_label"]
    ann.to_csv(os.path.join(OUT, "cluster_annotation_auto.tsv"), sep="\t")
    adata.obs["auto_label"] = adata.obs[CL].map(ann["auto_label"]).astype("category")
    log(f"自动注释：{adata.obs['auto_label'].value_counts().to_dict()}")
    rng = np.random.default_rng(SEED)
    idx = np.concatenate([rng.choice(v, min(len(v), 1000), replace=False) for v in adata.obs.groupby(CL, observed=True).indices.values()])
    sub = adata[np.sort(idx)].copy()
    sc.tl.rank_genes_groups(sub, CL, method="wilcoxon", n_genes=15)
    top = pd.DataFrame({c: list(sub.uns["rank_genes_groups"]["names"][c]) for c in sub.uns["rank_genes_groups"]["names"].dtype.names}).T
    top.index.name = CL
    top.to_csv(os.path.join(OUT, "cluster_top_markers.tsv"), sep="\t")
    for colr, fn in [(CL, "umap_leiden.png"), ("auto_label", "umap_auto_label.png"), (["tissue", "histology"], "umap_tissue_histology.png")]:
        fig = sc.pl.umap(adata, color=colr, show=False, return_fig=True, frameon=False,
                         legend_loc="on data" if colr == CL else "right margin")
        fig.savefig(os.path.join(OUT, fn), dpi=120, bbox_inches="tight"); plt.close(fig)
    adata.write_h5ad(os.path.join(OUT, f"{DS}_annotated.h5ad"), compression="lzf")
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
