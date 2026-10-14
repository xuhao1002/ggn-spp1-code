#!/usr/bin/env python
"""03 预备：列出配对表（不做任何差异分析）。AT2 谱系 = AT2 + AT2/Club-like epithelial + Proliferating epithelial。
前体 = 同患者所有 AAH/AIS/MIA 样本合并；LUAD = 同患者所有 LUAD 样本合并。
输出 results/03_pseudobulk/pairing_table.tsv、pairing_excluded.tsv、pairing_composition.tsv
用法：/opt/anaconda3/envs/ggn_target/bin/python scripts/03a_pairing_table.py"""
import os, numpy as np, pandas as pd, scanpy as sc
B = os.path.expanduser("~/全部数据库/GGN_target"); O = f"{B}/results/03_pseudobulk"
LIN = ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"]
a = sc.read_h5ad(f"{B}/results/01_snRNA/epithelial_cells.h5ad")
a = a[a.obs["final_label"].isin(LIN).values]
tot = np.asarray(a.layers["counts"].sum(axis=1)).ravel()
o = a.obs[["sample", "patient", "stage", "final_label"]].copy(); o["tot"] = tot
o["grp"] = o["stage"].map(lambda s: "LUAD" if s == "LUAD" else ("Precursor" if s in ("AAH", "AIS", "MIA") else "Normal"))
o = o[o["grp"] != "Normal"]
allsamp = pd.read_csv(f"{B}/results/01_snRNA/sample_meta.tsv", sep="\t")
rows, comp, excl = [], [], []
for p in sorted(allsamp["patient"].unique(), key=lambda x: int(x[1:])):
    sm = allsamp[allsamp["patient"] == p]
    prec = sm[sm["stage"].isin(["AAH", "AIS", "MIA"])]["sample"].tolist()
    prec_s = [f"{s}({sm.loc[sm['sample']==s,'stage'].iloc[0]})" for s in prec]
    lu = sm[sm["stage"] == "LUAD"]["sample"].tolist()
    d = o[o["patient"] == p]
    P, L = d[d["grp"] == "Precursor"], d[d["grp"] == "LUAD"]
    if not prec or not lu:
        excl.append(dict(patient=p, reason="无前体病变样本" if not prec else "无 LUAD 样本",
                         precursor=";".join(prec_s), luad=";".join(lu))); continue
    r = dict(patient=p, precursor_samples=";".join(prec_s), luad_samples=";".join(lu),
             precursor_stages="+".join(sorted(set(sm[sm['sample'].isin(prec)]['stage']))),
             n_cells_precursor=len(P), n_cells_luad=len(L), total_counts_precursor=int(P["tot"].sum()),
             total_counts_luad=int(L["tot"].sum()))
    for g, dd in [("precursor", P), ("luad", L)]:
        f = dd["final_label"].value_counts(normalize=True)
        for l, nm in zip(LIN, ["AT2", "ATClub", "Prolif"]):
            r[f"{nm}_frac_{g}"] = round(float(f.get(l, 0)), 3)
    # 敏感性（前体仅 AIS/MIA）
    P2 = P[P["stage"].isin(["AIS", "MIA"])]
    r["n_cells_precursor_AISMIA_only"] = len(P2)
    if len(P) < 50 or len(L) < 50:
        excl.append(dict(patient=p, reason=f"伪批量<50 细胞（前体 {len(P)}，LUAD {len(L)}）", precursor=";".join(prec_s), luad=";".join(lu)))
    rows.append(r)
t = pd.DataFrame(rows); t.to_csv(f"{O}/pairing_table.tsv", sep="\t", index=False)
pd.DataFrame(excl).to_csv(f"{O}/pairing_excluded.tsv", sep="\t", index=False)
c = pd.DataFrame({"AT2_luad_minus_prec": (t["AT2_frac_luad"] - t["AT2_frac_precursor"]).values,
                  "ATClub_luad_minus_prec": (t["ATClub_frac_luad"] - t["ATClub_frac_precursor"]).values,
                  "Prolif_luad_minus_prec": (t["Prolif_frac_luad"] - t["Prolif_frac_precursor"]).values}, index=t["patient"].values)
c.index.name = "patient"
c.round(3).to_csv(f"{O}/pairing_composition.tsv", sep="\t")
pd.set_option("display.width", 250, "display.max_columns", 30)
print(t.drop(columns=["precursor_samples", "luad_samples"]).to_string(index=False))
print(t[["patient", "precursor_samples", "luad_samples"]].to_string(index=False))
print("排除/提示：", pd.DataFrame(excl).to_string(index=False) if excl else "无")
print("配对数", len(t), "；仅 AIS/MIA 敏感性可用配对：",
      int(((t["n_cells_precursor_AISMIA_only"] >= 50)).sum()))
print("AT2 谱系三类占比（LUAD−前体）均值：", c.mean().round(3).to_dict())
print("LUAD 一侧：AT2/Club-like 占比中位", t["ATClub_frac_luad"].median(), "；前体一侧", t["ATClub_frac_precursor"].median())
