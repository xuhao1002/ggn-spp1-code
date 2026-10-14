#!/usr/bin/env python
"""生成第 3–5 层的候选表：并集 683 = 预设版 352（全部上皮验证）∪ AT2 谱系验证版 661。
set：both（330）/ only_main352（22）/ only_at2lineage（331）；direction = A 的 LUAD 对前体方向（上调=阻断靶点候选，下调=缺失/恢复清单）。
输出 results/07_candidates/candidates_683.tsv"""
import os, pandas as pd
B = os.path.expanduser("~/全部数据库/GGN_target/results")
a = pd.read_csv(f"{B}/06_validation/layer1_candidates.tsv", sep="\t", index_col=0)
b = pd.read_csv(f"{B}/06_validation_at2lineage/layer1_candidates.tsv", sep="\t", index_col=0)
u = pd.concat([a, b[~b.index.isin(a.index)]])
cols = ["log2FoldChange", "padj", "tier", "GSE164789_pass", "GSE189357_pass"]
c = u[cols].copy()
c["in_main352"] = c.index.isin(a.index); c["in_at2lineage661"] = c.index.isin(b.index)
c["set"] = c.apply(lambda r: "both" if r.in_main352 and r.in_at2lineage661 else ("only_main352" if r.in_main352 else "only_at2lineage"), axis=1)
c["direction"] = c["log2FoldChange"].map(lambda x: "up" if x > 0 else "down")
c.index.name = "gene"
c.to_csv(f"{B}/07_candidates/candidates_683.tsv", sep="\t")
print(len(c), c["set"].value_counts().to_dict(), c["direction"].value_counts().to_dict())
print("主 352 中：", c[c.in_main352]["direction"].value_counts().to_dict())
