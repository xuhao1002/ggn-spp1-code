"""A1d（方案 v1.3 §2-A1d）Visium 对照巨噬基因标记集敏感性。环境 ggn_spatial（与原 80_spatial_gene_exploratory.py 相同）。
样本读取、spot 过滤、归一化、平滑与纳入规则逐项沿用 80（GSE307534 部分）：组织 spot、filter_cells(min_counts=1)、
normalize_total 1e4 + log1p、spots ≥500 且平滑后 SPP1 非零 spot ≥100；M = 样本 SPP1 均值 ÷ 对照基因均值的均值；
同患者前体（AAH/AIS/MIA 取均值）vs 浸润（取均值）配对。
先用原 6 基因复现 E2（results/30_spatial/exploratory/E2_patient_pairs.tsv，21/24）；复现失败则不解读替代标记集。
替代标记集：去 CD68；去 C1QA/C1QB/C1QC；仅 MS4A6A+CD14。报告 "浸润 M > 前体 M" 的患者数相对 21/24 的变化，
并给出与 E2 相同的单侧 Wilcoxon 以便对照。
六个对照基因各自的患者级变化：按执行补充记录 §2（徐老师审定）——双侧 Wilcoxon 符号秩检验，六基因内 Holm（m=6），
作为探索性诊断报告；用于评估标志物的阶段相关变化，不作为筛选或替换对照基因的依据。
不显著不能证明标志物稳定；显著提示存在阶段相关变化，但不能单独判断来自细胞数量还是细胞内表达。
三套替代标记集全部报告，不得据这些 P 值择优。"""
import os, sys, glob, re, json, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump, holm_fixed
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy import stats
warnings.filterwarnings("ignore")

sha = check_plan(); log("方案哈希校验通过", sha[:16])

# 执行补充记录（A1c/A1d 实现细节，2026-10-19 执行前写定，不修改 v1.3）：哈希不符即拒绝运行
import hashlib as _hl
_SUPP = f"{P}/docs/执行补充_A1c_A1d_2026-10-19.md"
_SUPP_SHA = "e108a3449e1173b67292c7bf72a88d2b1a96af89c64ca8e8488e56a5c14707de"
if _hl.sha256(open(_SUPP, "rb").read()).hexdigest() != _SUPP_SHA:
    sys.exit("执行补充记录哈希不符，拒绝运行：" + _SUPP)
log("执行补充记录哈希校验通过", _SUPP_SHA[:16])
OUT = outdir("review_A1d")
CTRL = ["CD68", "C1QA", "C1QB", "C1QC", "MS4A6A", "CD14"]
SETS = {"original_6": CTRL,
        "drop_CD68": ["C1QA", "C1QB", "C1QC", "MS4A6A", "CD14"],
        "drop_C1Q": ["CD68", "MS4A6A", "CD14"],
        "MS4A6A_CD14_only": ["MS4A6A", "CD14"]}
OFFS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]

def adjacency(rc):                                   # 与 80 相同
    idx = {(r, c): i for i, (r, c) in enumerate(rc)}
    rows = []; cols = []
    for i, (r, c) in enumerate(rc):
        rows.append(i); cols.append(i)
        for dr, dc in OFFS:
            j = idx.get((r + dr, c + dc))
            if j is not None: rows.append(i); cols.append(j)
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(rc), len(rc)))
    d = np.asarray(A.sum(1)).ravel()
    return sp.diags(1 / d) @ A

rows = []
for d in sorted(glob.glob(f"{P}/raw/ST_GSE307534/extracted/GSM*")):
    base = os.path.basename(d); m = re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)(-\d+)?$", base)
    pat, stage = m.group(1), m.group(2)
    sub = glob.glob(d + "/*/filtered_feature_bc_matrix")[0]
    ad_ = sc.read_10x_mtx(sub, var_names="gene_symbols", make_unique=True)
    pos = pd.read_csv(os.path.dirname(sub) + "/spatial/tissue_positions.csv").set_index("barcode")
    pos = pos.loc[pos.index.intersection(ad_.obs_names)]; pos = pos[pos.in_tissue == 1]
    ad_ = ad_[pos.index].copy(); ad_.obsm["rc"] = pos[["array_row", "array_col"]].values.astype(int)
    sc.pp.filter_cells(ad_, min_counts=1)
    r = {"sample": base, "patient": pat, "stage": stage, "n_spots": int(ad_.n_obs)}
    if ad_.n_obs < 500:
        r["eligible"] = False; rows.append(r); continue
    sc.pp.normalize_total(ad_, target_sum=1e4); sc.pp.log1p(ad_)
    X = lambda g: np.asarray(ad_[:, g].X.todense()).ravel() if sp.issparse(ad_.X) else np.asarray(ad_[:, g].X).ravel()
    spp1 = X("SPP1"); W = adjacency(ad_.obsm["rc"])
    r["spp1_nonzero_smoothed"] = int(((W @ spp1) > 0).sum())
    r["eligible"] = bool(r["spp1_nonzero_smoothed"] >= 100)
    r["SPP1_mean"] = float(spp1.mean())
    for g in CTRL: r[f"mean_{g}"] = float(X(g).mean()) if g in ad_.var_names else np.nan
    rows.append(r); log(base, r["n_spots"], r["eligible"])
S = pd.DataFrame(rows)
S["group"] = S.stage.map({"AAH": "precursor", "AIS": "precursor", "MIA": "precursor", "LUAD": "LUAD", "Normal": "Normal"})
S.to_csv(OUT + "A1d_sample_level.tsv", sep="\t", index=False)
E = S[S.eligible.fillna(False).astype(bool)].copy()
for name, gs in SETS.items():
    cm = E[[f"mean_{g}" for g in gs]].mean(axis=1)
    E[f"M_{name}"] = np.where(cm > 0, E.SPP1_mean / cm, np.nan)
cols = [c for c in E.columns if c.startswith("M_") or c.startswith("mean_")] + ["SPP1_mean"]
pre = E[E.group == "precursor"].groupby("patient")[cols].mean(); lu = E[E.group == "LUAD"].groupby("patient")[cols].mean()
pairs = pre.join(lu, lsuffix="_pre", rsuffix="_lu", how="inner")
pairs.to_csv(OUT + "A1d_patient_pairs.tsv", sep="\t")

# ---- 复现核对（与原 E2 相同的 M） ----
orig = pd.read_csv(f"{P}/results/30_spatial/exploratory/E2_patient_pairs.tsv", sep="\t", index_col=0)
j = orig.join(pairs[["M_original_6_pre", "M_original_6_lu"]], how="inner")
repro = bool(len(j) == len(orig) == len(pairs)
             and np.allclose(j.M_pre, j.M_original_6_pre, rtol=1e-6, atol=1e-9)
             and np.allclose(j.M_lu, j.M_original_6_lu, rtol=1e-6, atol=1e-9))
report = {"plan_sha256": sha, "reproduces_original_E2": repro, "n_pairs": int(len(pairs)),
          "original_reported": "E2: 21/24 patients higher in invasive; one-sided Wilcoxon P 1.5e-5 (Holm 4.5e-5)"}
res = {}
for name in SETS:
    d = (pairs[f"M_{name}_lu"] - pairs[f"M_{name}_pre"]).dropna().values
    t = stats.wilcoxon(d, alternative="greater") if len(d) and not np.all(d == 0) else None
    res[name] = {"genes": SETS[name], "n_pairs": int(len(d)), "n_LUAD_higher": int((d > 0).sum()),
                 "median_diff": float(np.median(d)) if len(d) else np.nan,
                 "wilcoxon_one_sided_p_same_test_as_E2": float(t.pvalue) if t is not None else None}
report["marker_sets"] = res
# ---- 六个对照基因各自的患者级变化 ----
pv, genes = {}, {}
for g in CTRL:
    d = (pairs[f"mean_{g}_lu"] - pairs[f"mean_{g}_pre"]).dropna().values
    t = stats.wilcoxon(d, alternative="two-sided") if len(d) and not np.all(d == 0) else None
    genes[g] = {"n_pairs": int(len(d)), "n_LUAD_higher": int((d > 0).sum()), "median_diff": float(np.median(d)) if len(d) else np.nan,
                "mean_diff": float(np.mean(d)) if len(d) else np.nan,
                "wilcoxon_two_sided_p": float(t.pvalue) if t is not None else None}
    pv[g] = genes[g]["wilcoxon_two_sided_p"]
for g, v in holm_fixed(pv, len(CTRL)).items(): genes[g]["holm_m6"] = v
for g in CTRL: genes[g]["direction"] = ("higher in invasive" if genes[g]["median_diff"] > 0 else ("lower in invasive" if genes[g]["median_diff"] < 0 else "no change"))
report["control_genes_per_gene"] = genes
report["control_genes_rule"] = ("执行补充记录 §2：六个对照基因的患者级配对变化采用双侧 Wilcoxon 符号秩检验，并在六个基因内进行 Holm 校正；"
                                "报告效应方向、变化量、有效配对数、原始及校正后 P；用于评估标志物的阶段相关变化，不作为筛选或替换对照基因的依据。"
                                "不显著不能证明标志物稳定；显著提示存在阶段相关变化，但不能单独判断变化来自细胞数量还是细胞内表达。"
                                "三套替代标记集全部报告，不得据这些 P 值择优。若对照基因本身随分期系统变化，须在稿件中写明（方案 A1d）。")
report["supplementary_record_sha256"] = _SUPP_SHA
report["note"] = ("复现失败：先查原因，不解读替代标记集。" if not repro else "复现成功：替代标记集结果可解读。")
dump(report, OUT + "A1d_results.json")
log("A1d 完成", {k: (v["n_LUAD_higher"], v["n_pairs"]) for k, v in res.items()})
