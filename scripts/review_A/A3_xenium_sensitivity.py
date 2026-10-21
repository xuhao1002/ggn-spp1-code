"""A3（方案 v1.3 §2-A3、§6）Xenium 来源归属敏感性。环境 ggn_target。
无转录本坐标/边界/图像，故不做重新分割；仅做替代分型规则等敏感性。
所有来源比例一律限定为 "under the specified spillover model"。"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy import stats
from scipy.spatial import cKDTree

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A3")
RAW = f"{P}/raw/Xenium_Takano2024"

# 与 103_xenium_spillover.py 完全一致的标记集与参考基因
# 标记集逐字取自 103_xenium_spillover.py 第 8-9 行，不另行定义，以保证与原分析一致。
SETS = {"Macrophage": ["CD68","CD163","MARCO","MSR1","C1QC","MPEG1","MS4A4A"],
        "Epithelial": ["KRT8","NAPSA","SFTPC","SFTPB","NKX2-1","CDH1","SCGB3A2"],
        "Fibroblast": ["COL1A1","PDGFRA","FAP","POSTN"],
        "Endothelial": ["PECAM1","KDR","PLVAP"],
        "T/NK": ["CD3E","CD8A","GZMA","LCK"],
        "B/Plasma": ["CD19","MZB1","BANK1"]}
# 断言：与 103 的定义逐字一致，否则中止（防止两处漂移）
import re as _re
_src = open(f"{P}/scripts/103_xenium_spillover.py", encoding="utf-8").read()
for _k, _v in SETS.items():
    for _g in _v:
        assert f'"{_g}"' in _src, f"标记基因 {_g} 未出现在 103_xenium_spillover.py 中，分型规则可能已漂移"
REF = ["CD68", "C1QC", "MARCO", "MSR1"]
EPI_REV = ["EPCAM", "KRT8", "SFTPC"]          # A3c 反向标记混入指标（诊断性，非"反向泄漏率"）
RULES = {"base_ge2": ("thresh", 2), "thresh_ge3": ("thresh", 3),
         "thresh_ge5": ("thresh", 5), "ratio_2x": ("ratio", 2.0)}

def assign(M, names, rule):
    o = np.argsort(-M, axis=1)
    top = M[np.arange(len(M)), o[:, 0]]; sec = M[np.arange(len(M)), o[:, 1]]
    if rule[0] == "thresh":
        ok = (top >= rule[1]) & (top > sec)
    else:
        ok = (top >= 2) & (top >= rule[1] * np.maximum(sec, 1e-9))
    return np.where(ok, names[o[:, 0]], "Unassigned")

rows = []
import glob, re as _re2
for dd in sorted(glob.glob(RAW + "/output-*")):
    sid = _re2.search(r"__(LUAD_No\d+|TSU-\d+)__", dd).group(1)
    grp = "invasive" if sid.startswith("LUAD") else "AIS/MIA"
    a = sc.read_10x_h5(f"{dd}/cell_feature_matrix.h5"); a.var_names_make_unique()
    cells = pd.read_csv(f"{dd}/cells.csv.gz").set_index("cell_id"); cells.index = cells.index.astype(str)
    a = a[a.obs_names.isin(cells.index)].copy()
    xy = cells.loc[a.obs_names, ["x_centroid", "y_centroid"]].values
    X = sp.csr_matrix(a.X); gi = lambda g: a.var_names.get_loc(g)
    M = np.column_stack([np.asarray(X[:, [gi(g) for g in v if g in a.var_names]].sum(1)).ravel()
                         for v in SETS.values()]); names = np.array(list(SETS))
    spp1 = np.asarray(X[:, gi("SPP1")].todense()).ravel()
    ref = np.asarray(X[:, [gi(g) for g in REF if g in a.var_names]].sum(1)).ravel()
    rev = np.asarray(X[:, [gi(g) for g in EPI_REV if g in a.var_names]].sum(1)).ravel()
    tot = spp1.sum()

    for rname, rule in RULES.items():
        ct = assign(M, names, rule)
        mac = ct == "Macrophage"; epi = ct == "Epithelial"; una = ct == "Unassigned"
        r = {"section": sid, "group": grp, "rule": rname, "ref_set": "all4",
             "n_cells": int(len(ct)), "n_mac": int(mac.sum()), "n_epi": int(epi.sum()),
             # A3a：Unassigned 的细胞比例，以及其承载的 SPP1 转录本比例（后者此前未报告）
             "unassigned_cell_frac": float(una.mean()),
             "unassigned_SPP1_frac": float(spp1[una].sum() / tot) if tot > 0 else np.nan}
        if mac.sum() == 0 or tot == 0:
            r.update(leak_L=np.nan, obs_over_exp=np.nan, epi_share_raw=np.nan,
                     epi_share_corrected=np.nan, mac_share_raw=np.nan, mac_share_corrected=np.nan,
                     evaluable=False); rows.append(r); continue
        L = ref[epi].sum() / ref[mac].sum(); exp_ = L * spp1[mac].sum(); obs = spp1[epi].sum()
        move = min(exp_, obs)
        r.update(leak_L=float(L), obs_over_exp=float(obs / exp_) if exp_ > 0 else np.nan,
                 epi_share_raw=float(obs / tot), epi_share_corrected=float((obs - move) / tot),
                 mac_share_raw=float(spp1[mac].sum() / tot),
                 mac_share_corrected=float((spp1[mac].sum() + move) / tot), evaluable=True)
        # A3c 反向标记混入指标（诊断性：巨噬细胞内上皮特征基因的相对计数）
        r["reverse_admixture_index"] = float(rev[mac].sum() / rev[epi].sum()) if rev[epi].sum() > 0 else np.nan
        # locality：仅在 SPP1+ 巨噬 >=20 的切片执行，定义同 103（上皮自身 SPP1 vs 到最近 SPP1+ 巨噬距离）
        pos = mac & (spp1 >= 2)
        if pos.sum() >= 20 and epi.sum() >= 10:
            dist = cKDTree(xy[pos]).query(xy[epi])[0]
            rho, pv = stats.spearmanr(spp1[epi], dist)
            r["locality_rho"] = float(rho); r["locality_p"] = float(pv)
        rows.append(r)

    # A3b 参考基因留一 + 仅 CD68（固定在 base 规则下）
    ct = assign(M, names, RULES["base_ge2"]); mac = ct == "Macrophage"; epi = ct == "Epithelial"
    for sub_name, sub in [(f"drop_{g}", [x for x in REF if x != g]) for g in REF] + [("CD68_only", ["CD68"])]:
        rr = np.asarray(X[:, [gi(g) for g in sub if g in a.var_names]].sum(1)).ravel()
        if mac.sum() == 0 or rr[mac].sum() == 0 or tot == 0: continue
        L = rr[epi].sum() / rr[mac].sum(); exp_ = L * spp1[mac].sum(); obs = spp1[epi].sum()
        move = min(exp_, obs)
        rows.append({"section": sid, "group": grp, "rule": "base_ge2", "ref_set": sub_name,
                     "n_cells": int(len(ct)), "n_mac": int(mac.sum()), "n_epi": int(epi.sum()),
                     "leak_L": float(L), "obs_over_exp": float(obs / exp_) if exp_ > 0 else np.nan,
                     "epi_share_raw": float(obs / tot), "epi_share_corrected": float((obs - move) / tot),
                     "mac_share_raw": float(spp1[mac].sum() / tot),
                     "mac_share_corrected": float((spp1[mac].sum() + move) / tot), "evaluable": True})
    log("完成切片", sid, grp)

T = pd.DataFrame(rows); T.to_csv(OUT + "A3_per_section.tsv", sep="\t", index=False)

# 决策规则：任一替代设定下，浸润组校正后巨噬来源比例中位数跌破 50% → 删去 "main source"
summ = {}
for (rule, refset), g in T[T.evaluable].groupby(["rule", "ref_set"]):
    for grp, gg in g.groupby("group"):
        summ[f"{rule}|{refset}|{grp}"] = {
            "n_sections": int(len(gg)),
            "median_mac_share_corrected": float(gg.mac_share_corrected.median()),
            "median_epi_share_corrected": float(gg.epi_share_corrected.median()),
            "median_obs_over_exp": float(gg.obs_over_exp.median()),
            "median_unassigned_cell_frac": float(gg.unassigned_cell_frac.median()) if "unassigned_cell_frac" in gg else None,
            "median_unassigned_SPP1_frac": float(gg.unassigned_SPP1_frac.median()) if "unassigned_SPP1_frac" in gg else None}
inv = {k: v for k, v in summ.items() if k.endswith("|invasive")}
below = {k: v["median_mac_share_corrected"] for k, v in inv.items() if v["median_mac_share_corrected"] < 0.5}
dump({"plan_sha256": sha, "summary": summ,
      "decision_rule": "任一替代设定下浸润组校正后巨噬来源中位数 <50% → 正文删去 'macrophages were its main source'，改报区间",
      "settings_below_50pct": below,
      "decision": ("删去 main source，改报区间" if below else "保留 main source 的表述"),
      "note": ("反向不成立：所有设定都高于 50% 不构成巨噬细胞为真实主要来源的证据，"
               "来源比例一律限定于本模型（方案 A3.0）。")}, OUT + "A3_results.json")
log("A3 完成")
