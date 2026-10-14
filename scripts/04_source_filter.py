#!/usr/bin/env python
"""
来源过滤（改规则版，2026-10-08 用户确认；旧版见 04_source_filter_v1.py）：
- 谱系：上皮 + 8 个非上皮谱系：T/NK、B/Plasma、Myeloid、Mast、Fibroblast、SMC/pericyte、Endothelial、Lymphatic EC
- 每个谱系、每个样本：求平均 log1p 表达（该谱系在该样本 >= 20 个细胞才计入），再取样本中位数
- 比值 = 上皮 / 表达最高的非上皮谱系；比值 < 0.2 → 标记
用法：
    python scripts/04_source_filter.py peng        # 用 snRNA_final.h5ad（final_label）
    python scripts/04_source_filter.py GSE164789   # 用 05 预处理并注释后的对象（cluster_annotation_auto.tsv 的 final_label）
    python scripts/04_source_filter.py rescue      # 汇总 Peng 与 GSE164789：被救回基因、A 前 50 被标记基因及来源谱系
输出 results/04_source_filter/：ratio_<数据集>.tsv、flagged_<数据集>.txt、rescued_genes.tsv、A_top50_source.tsv、source_filter_final.tsv
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
B = os.path.expanduser("~/全部数据库/GGN_target"); O = f"{B}/results/04_source_filter"
os.makedirs(O, exist_ok=True)
MIN_CELLS = 20
LINEAGE_OF = {"T cell": "T/NK", "NK": "T/NK", "Proliferating T/NK": "T/NK",
              "B cell": "B/Plasma", "Plasma": "B/Plasma", "Proliferating B/Plasma": "B/Plasma",
              "Macrophage": "Myeloid", "Monocyte": "Myeloid", "cDC": "Myeloid", "pDC": "Myeloid",
              "Neutrophil": "Myeloid", "Proliferating myeloid": "Myeloid",
              "Mast": "Mast", "Fibroblast": "Fibroblast", "SMC/pericyte": "SMC/pericyte",
              "Endothelial": "Endothelial", "Lymphatic EC": "Lymphatic EC"}
EPI = ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal", "Proliferating epithelial",
       "AT2/Club-like epithelial", "AT2-like (patient-specific)"]

def ratio_table(a, label_col, sample_col):
    lab = a.obs[label_col].astype(str).values
    smp = a.obs[sample_col].astype(str).values
    lin = np.array([("Epithelial" if l in EPI else LINEAGE_OF.get(l, "")) for l in lab])
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    med, ncell = {}, {}
    for L in ["Epithelial"] + sorted(set(LINEAGE_OF.values())):
        rows = []
        for s in np.unique(smp):
            idx = np.where((lin == L) & (smp == s))[0]
            if len(idx) >= MIN_CELLS:
                rows.append(np.asarray(X[idx].mean(axis=0)).ravel())
        ncell[L] = len(rows)
        if rows:
            med[L] = np.median(np.vstack(rows), axis=0)
    M = pd.DataFrame(med, index=a.var_names)
    non = [c for c in M.columns if c != "Epithelial"]
    mx = M[non].max(axis=1); arg = M[non].idxmax(axis=1)
    r = pd.DataFrame({"epi": M["Epithelial"], "max_nonepi": mx, "max_nonepi_lineage": arg})
    r["ratio"] = r["epi"] / r["max_nonepi"].replace(0, np.nan)
    r.loc[(r["max_nonepi"] == 0) & (r["epi"] > 0), "ratio"] = np.inf
    r["flag"] = r["ratio"] < 0.2          # 都不表达（比值 NaN）不判定
    r = r.join(M[non].add_prefix("med_"))
    r.index.name = "gene"
    return r, ncell

def main():
    mode = sys.argv[1]
    if mode == "peng":
        a = sc.read_h5ad(f"{B}/results/01_snRNA/snRNA_final.h5ad")
        r, nc = ratio_table(a, "final_label", "sample")
        r.to_csv(f"{O}/ratio_peng.tsv", sep="\t")
        open(f"{O}/flagged_peng.txt", "w").write("\n".join(r.index[r["flag"]]))
        print("各谱系有效样本数：", nc, "；被标记基因", int(r["flag"].sum()), "/", len(r))
    elif mode == "GSE164789":
        d = f"{B}/results/05_validation/GSE164789"
        a = sc.read_h5ad(f"{d}/GSE164789_annotated.h5ad")
        ann = pd.read_csv(f"{d}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
        a.obs["final_label"] = a.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
        r, nc = ratio_table(a, "final_label", "sample")
        r.to_csv(f"{O}/ratio_GSE164789.tsv", sep="\t")
        open(f"{O}/flagged_GSE164789.txt", "w").write("\n".join(r.index[r["flag"]]))
        print("各谱系有效样本数：", nc, "；被标记基因", int(r["flag"].sum()), "/", len(r))
    else:
        P = pd.read_csv(f"{O}/ratio_peng.tsv", sep="\t", index_col=0)
        G = pd.read_csv(f"{O}/ratio_GSE164789.tsv", sep="\t", index_col=0)
        j = P[["ratio", "max_nonepi_lineage", "flag"]].join(G[["ratio", "max_nonepi_lineage"]], rsuffix="_GSE164789", how="left")
        j = j.rename(columns={"flag": "flag_peng", "ratio": "ratio_peng", "max_nonepi_lineage": "lineage_peng",
                              "max_nonepi_lineage_GSE164789": "lineage_GSE164789"})
        j["rescued"] = j["flag_peng"] & (j["ratio_GSE164789"] >= 0.2)
        j["final_flag"] = j["flag_peng"] & ~j["rescued"]
        j.index.name = "gene"
        j.to_csv(f"{O}/source_filter_final.tsv", sep="\t")
        j[j["rescued"]].to_csv(f"{O}/rescued_genes.tsv", sep="\t")
        A = pd.read_csv(f"{B}/results/03_pseudobulk/deg_A_main.tsv", sep="\t", index_col=0)
        sig = A[A["padj"] < 0.05].join(j, how="left")
        sig.sort_values("padj").head(50).to_csv(f"{O}/A_top50_source.tsv", sep="\t")
        print("Peng 被标记", int(j["flag_peng"].sum()), "；被救回", int(j["rescued"].sum()), "；最终标记", int(j["final_flag"].sum()))
        print("A 显著基因", len(sig), "；Peng 标记", int(sig["flag_peng"].sum()), "；最终标记", int(sig["final_flag"].sum()),
              "；被救回", int(sig["rescued"].sum()), "；GSE164789 中无此基因", int(sig["ratio_GSE164789"].isna().sum()))
        t = sig.sort_values("padj").head(50)
        print("A 前 50 被标记（Peng）：")
        for g, x in t[t["flag_peng"]].iterrows():
            print(f"  {g}: Peng ratio={x.ratio_peng:.2f}（{x.lineage_peng}），GSE164789 ratio={x.ratio_GSE164789:.2f}（{x.lineage_GSE164789}），"
                  f"{'被救回' if x.rescued else '仍标记'}，LFC{x.log2FoldChange:+.1f}")
        print("关注基因：")
        for g in ["CTHRC1", "MMP11", "VCAN", "SPP1", "CXCL14", "ARHGEF1", "IGFBP3", "FAM83A", "SLC2A1", "PLAU", "CORO1A", "LCP1", "CD53"]:
            if g in j.index:
                x = j.loc[g]; print(f"  {g}: Peng {x.ratio_peng:.2f}（{x.lineage_peng}）{'标记' if x.flag_peng else ''}；GSE164789 {x.ratio_GSE164789:.2f}；{'被救回' if x.rescued else ''}")

if __name__ == "__main__":
    main()
