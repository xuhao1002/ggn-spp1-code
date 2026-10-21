"""A2b（方案 v1.3 §2-A2b）：两个外部单细胞队列的髓系 SPP1 单位级散点（按组分列，叠加中位数）。描述性，无检验（A2c）。
输入：results/review_A2/A2_units_GSE164789.tsv、A2_units_GSE189357.tsv（A2 输出，未改动）。
GSE164789 中患者 HYY 贡献两个单位（偏离记录 D2）：两点以细线相连并标注。"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A2b")
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42, "savefig.dpi": 300})
PRE, INV, INK, MUT = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
rng = np.random.default_rng(20261019)
fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.6), gridspec_kw={"wspace": 0.45})
rows = []
for ax, ds in zip(axes, ["GSE164789", "GSE189357"]):
    u = pd.read_csv(f"{P}/results/review_A2/A2_units_{ds}.tsv", sep="\t")
    u = u[u.included == True].copy()
    u["x"] = u.condition.map({"Precursor": 0, "IAC": 1}); assert u.x.notna().all()
    u["xj"] = u.x + rng.uniform(-0.13, 0.13, len(u))
    for x, c in [(0, PRE), (1, INV)]:
        v = u[u.x == x]
        ax.scatter(v.xj, v.SPP1_pseudobulk_log1p_cp10k, s=14, color=c, edgecolor="white", lw=0.4, zorder=3)
        med = v.SPP1_pseudobulk_log1p_cp10k.median(); ax.plot([x - 0.25, x + 0.25], [med, med], color=INK, lw=1.3, zorder=4)
        rows.append(dict(cohort=ds, group=["Precursor", "IAC"][x], n_units=len(v), median=med))
    multi = u.groupby("patient").filter(lambda g: g.condition.nunique() > 1)
    for p, g in multi.groupby("patient"):
        g = g.sort_values("x"); ax.plot(g.xj, g.SPP1_pseudobulk_log1p_cp10k, color=MUT, lw=0.7, ls=(0, (2, 1.5)), zorder=2)
        r = g.iloc[-1]; ax.annotate(f"{p} (2 units)", (r.xj, r.SPP1_pseudobulk_log1p_cp10k), xytext=(6, 0), textcoords="offset points",
                                    fontsize=6, color=MUT, va="center")
    n0, n1 = (u.x == 0).sum(), (u.x == 1).sum()
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"Precursor\n(n = {n0})", f"IAC\n(n = {n1})"]); ax.set_xlim(-0.5, 1.6)
    ax.set_ylabel("Myeloid SPP1, log1p(CP10k)"); ax.set_title(f"{ds} ({u.patient.nunique()} patients)", loc="left", fontsize=7.5)
    u.drop(columns=["xj"]).to_csv(OUT + f"A2b_plotted_units_{ds}.tsv", sep="\t", index=False)
pd.DataFrame(rows).to_csv(OUT + "A2b_group_medians.tsv", sep="\t", index=False)
for e in ["png", "pdf", "tiff"]: fig.savefig(OUT + f"FigureS1_A2b_unit_scatter.{e}", bbox_inches="tight", dpi=300)
log("A2b 完成", rows)
