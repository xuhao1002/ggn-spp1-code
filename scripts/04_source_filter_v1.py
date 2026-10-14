#!/usr/bin/env python
"""来源过滤：每个基因 上皮平均表达 / 表达最高的非上皮类型平均表达（snRNA_final，log1p 归一化 X）。
上皮 = final_compartment == Epithelial；非上皮类型 = final_label 细胞数 >= 500 的非上皮类；比值 < 0.2 → suspected_nonepithelial。
输出 results/04_source_filter/gene_source_ratio.tsv、flagged_genes.txt、A_top50_source.tsv
用法：/opt/anaconda3/envs/ggn_target/bin/python scripts/04_source_filter.py"""
import os, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
B = os.path.expanduser("~/全部数据库/GGN_target"); O = f"{B}/results/04_source_filter"
a = sc.read_h5ad(f"{B}/results/01_snRNA/snRNA_final.h5ad")
lab = a.obs["final_label"].astype(str).values; comp = a.obs["final_compartment"].astype(str).values
epi = comp == "Epithelial"
types = [l for l in pd.unique(lab[~epi]) if (lab == l).sum() >= 500]
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
def mean_of(mask): return np.asarray(X[np.where(mask)[0]].mean(axis=0)).ravel()
m_epi = mean_of(epi)
M = pd.DataFrame({t: mean_of(lab == t) for t in types}, index=a.var_names)
mx = M.max(axis=1); arg = M.idxmax(axis=1)
r = pd.DataFrame({"epi_mean": m_epi, "max_nonepi_mean": mx.values, "max_nonepi_type": arg.values}, index=a.var_names)
r["ratio"] = r["epi_mean"] / r["max_nonepi_mean"].replace(0, np.nan)
r.loc[r["max_nonepi_mean"] == 0, "ratio"] = np.inf
r["suspected_nonepithelial"] = r["ratio"] < 0.2
r.index.name = "gene"
r.to_csv(f"{O}/gene_source_ratio.tsv", sep="\t")
open(f"{O}/flagged_genes.txt", "w").write("\n".join(r.index[r["suspected_nonepithelial"]]))
A = pd.read_csv(f"{B}/results/03_pseudobulk/deg_A_main.tsv", sep="\t", index_col=0)
sig = A[A["padj"] < 0.05]
j = sig.join(r, how="left")
j.sort_values("padj").head(50).to_csv(f"{O}/A_top50_source.tsv", sep="\t")
print("非上皮类型（>=500 细胞）：", types)
print("全部基因被标记：", int(r["suspected_nonepithelial"].sum()), "/", len(r))
print("A 显著基因被标记：", int(j["suspected_nonepithelial"].sum()), "/", len(j),
      "；其中上调", int((j["suspected_nonepithelial"] & (j["log2FoldChange"] > 0)).sum()),
      "下调", int((j["suspected_nonepithelial"] & (j["log2FoldChange"] < 0)).sum()))
t50 = j.sort_values("padj").head(50)
print("A 前 50 中被标记：", ", ".join(f"{g}(ratio={x.ratio:.2f},{x.max_nonepi_type},LFC{x.log2FoldChange:+.1f})" for g, x in t50[t50["suspected_nonepithelial"]].iterrows()) or "无")
for g in ["CTHRC1", "MMP11", "VCAN", "SPP1", "CXCL14", "ARHGEF1", "IGFBP3", "FAM83A", "SLC2A1", "PLAU", "ITGB8", "FHL2", "CORO1A", "LCP1", "CD53", "SFTPC", "SFTPB"]:
    if g in r.index: print(g, "ratio=%.2f" % r.loc[g, "ratio"], r.loc[g, "max_nonepi_type"], "flag" if r.loc[g, "suspected_nonepithelial"] else "")
