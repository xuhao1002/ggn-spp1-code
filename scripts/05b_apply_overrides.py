#!/usr/bin/env python
"""05 自动注释后的人工复核（2026-10-08，依据 cluster_annotation_auto.tsv 的 MKI67/PTPRC/EPCAM 均值与得分差）。
GSE189357：簇 11 → Proliferating myeloid（MKI67 0.855，得分 0.185）；簇 27 → Doublet/ambiguous（巨噬与内皮得分几乎相同，219 细胞/2 患者）
GSE164789：簇 7 → Proliferating myeloid（MKI67 1.048，得分 0.183）；簇 25 → Doublet/ambiguous（巨噬+AT2 共表达，109 细胞）
其余簇沿用自动标签。用法：python scripts/05b_apply_overrides.py"""
import os, pandas as pd
B = os.path.expanduser("~/全部数据库/GGN_target/results/05_validation")
OV = {"GSE189357": {11: "Proliferating myeloid", 27: "Doublet/ambiguous"},
      "GSE164789": {7: "Proliferating myeloid", 25: "Doublet/ambiguous"}}
for ds, ov in OV.items():
    p = f"{B}/{ds}/cluster_annotation_auto.tsv"
    a = pd.read_csv(p, sep="\t", index_col=0)
    a["final_label"] = a["auto_label"]
    for c, l in ov.items(): a.loc[c, "final_label"] = l
    a.to_csv(p, sep="\t")
    print(ds, ov)
