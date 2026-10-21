"""A4（方案 v1.3 §2-A4、§6.1、§6.2）分期合并的敏感性。环境 ggn_target。
A4.0（标题与时间性表述）已无条件落实于稿件，本脚本不涉及。
A4a  仅 AIS vs 浸润，髓系 SPP1 配对伪批量检验（与原主分析同一管线），<8 对记为不可评估。只用于说明方向与量级，不决定标题。
A4b  AAH 单独描述：逐患者列出数值，不做推断。
A4c  病灶等权 vs 细胞数加权（§6.2）：同一组合格病灶、同一组 x_j、同一患者集合、同一检验，只差权重。
     A1a（负二项 ~ patient + condition，聚合原始计数）为原始主模型，在 A1 脚本中单独报告，不在此与 A4c 比较。"""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump, boot_ci, SEED
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
warnings.filterwarnings("ignore")

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A4")
MIN_CELLS = 50; MIN_PAIRS = 8
PRE = ("AAH", "AIS", "MIA"); INV = "LUAD"
MYE = ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"]   # 与 14 的 Myeloid 一致
ORIG_MAIN = {"log2FC": 3.5747252005501227, "pvalue": 8.31924610549111e-12, "padj": 3.6584601925285563e-09}  # 14 的髓系 SPP1，供对照

# ---- 与 14_microenv_layer1.py 完全相同 ----
def filt(c):
    n = len(c); k = ((c >= 10).sum(axis=0) >= np.ceil(n / 3)) & (c.sum(axis=0) >= 50); return c.loc[:, k.values]

def deseq(counts, meta, design, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    c = filt(counts)
    dds = DeseqDataSet(counts=c, metadata=meta, design=design, refit_cooks=True, n_cpus=6, quiet=True); dds.deseq2()
    ds = DeseqStats(dds, contrast=contrast, quiet=True); ds.summary()
    r = ds.results_df.copy(); r.index.name = "gene"; return r

a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad"); log("载入", a.shape)
assert "counts" in a.layers, "缺少 layers['counts']（.X 为对数归一化值，不可用）"
C = a.layers["counts"].tocsr()
gi = {g: i for i, g in enumerate(a.var_names)}
ob = a.obs[["patient", "sample", "stage", "final_label"]].astype(str).copy()
genes = a.var_names.copy()
del a
ep = np.asarray(C[:, gi["EPCAM"]].todense()).ravel()
spp1_all = np.asarray(C[:, gi["SPP1"]].todense()).ravel()
lib_all = np.asarray(C.sum(axis=1)).ravel()                 # 全基因总计数（§6.1 分母）
report = {"plan_sha256": sha}

# ======================= A4a：AIS vs 浸润，髓系配对伪批量 =======================
mye = ob["final_label"].isin(MYE).values & (ep < 2)
groups, meta, excl = {}, [], []
for p in sorted(ob.patient.unique(), key=lambda x: int(x[1:])):
    idx = {}
    for g, stages in [("Precursor", ("AIS",)), ("LUAD", (INV,))]:
        idx[g] = np.where(mye & (ob.patient.values == p) & ob.stage.isin(stages).values)[0]
    if min(len(idx["Precursor"]), len(idx["LUAD"])) < MIN_CELLS:
        excl.append(p); continue
    for g in ["Precursor", "LUAD"]:
        groups[f"{p}_{g}"] = idx[g]; meta.append(dict(unit=f"{p}_{g}", patient=p, condition=g, n_cells=len(idx[g])))
meta = pd.DataFrame(meta).set_index("unit") if meta else pd.DataFrame()
n_pairs = len(meta) // 2 if len(meta) else 0
a4a = {"definition": "髓系（同 14）且 EPCAM<2；前体仅 AIS（排除 AAH 与 MIA）；每侧 ≥50 细胞；pydeseq2 ~ patient + condition；LUAD vs Precursor",
       "n_pairs": int(n_pairs), "patients": sorted(meta.patient.unique().tolist()) if n_pairs else [],
       "patients_not_pairable": excl, "reference_all_precursor": ORIG_MAIN}
if n_pairs < MIN_PAIRS:
    a4a.update(evaluable=False, reason=f"可配对患者 {n_pairs} < {MIN_PAIRS}")
    log("A4a 不可评估", n_pairs)
else:
    counts = pd.DataFrame({k: np.asarray(C[i].sum(axis=0)).ravel() for k, i in groups.items()}, index=genes).T.astype(int)
    m = meta.copy(); m["condition"] = pd.Categorical(m["condition"], categories=["Precursor", "LUAD"])
    r = deseq(counts, m[["patient", "condition"]], "~ patient + condition", ["condition", "LUAD", "Precursor"])
    s = r.loc["SPP1"]
    a4a.update(evaluable=True, SPP1_log2FC=float(s.log2FoldChange), SPP1_lfcSE=float(s.lfcSE),
               SPP1_p_two_sided=float(s.pvalue), SPP1_padj=float(s.padj),
               SPP1_p_one_sided=float(s.pvalue / 2 if s.log2FoldChange > 0 else 1 - s.pvalue / 2),
               direction=("same as main analysis" if s.log2FoldChange > 0 else "reversed"),
               use="仅说明分期合并是否影响方向与量级；不用于决定标题（方案 A4）")
    m.to_csv(OUT + "A4a_units.tsv", sep="\t")
    log("A4a", n_pairs, "对", {k: a4a[k] for k in ["SPP1_log2FC", "SPP1_padj"]})
report["A4a"] = a4a

# ======================= 严格巨噬细胞病灶表（A4b、A4c 共用 x_j） =======================
mac = (ob["final_label"].values == "Macrophage") & (ep < 2)          # 与 A1a 同一筛选
les = []
for (p, smp, stg), idx in ob[mac].groupby(["patient", "sample", "stage"]).indices.items():
    rows_ = np.where(mac)[0][idx]
    n = len(rows_); tot = lib_all[rows_].sum(); s1 = spp1_all[rows_].sum()
    les.append(dict(patient=p, lesion=smp, stage=stg,
                    group=("Precursor" if stg in PRE else ("LUAD" if stg == INV else "Normal")),
                    n_mac=int(n), SPP1_counts=int(s1), total_counts=int(tot),
                    x=float(np.log1p(s1 / tot * 1e4)) if tot > 0 else np.nan,     # §6.1：先聚合，再 CP10k，再 log1p
                    qualifies=bool(n >= MIN_CELLS)))
L = pd.DataFrame(les).sort_values(["patient", "group", "lesion"])
L.to_csv(OUT + "A4_lesions_strict_macrophage.tsv", sep="\t", index=False)

# ======================= A4b：AAH 逐患者描述（不推断） =======================
aah_pat = sorted(L[L.stage == "AAH"].patient.unique(), key=lambda x: int(x[1:]))
b = []
for p in aah_pat:
    sub = L[L.patient == p]
    for _, r in sub[sub.stage.isin(["AAH", INV, "Normal"])].iterrows():
        b.append(dict(patient=p, lesion=r.lesion, stage=r.stage, n_mac=r.n_mac, SPP1_log1p_cp10k=r.x, qualifies=r.qualifies))
B = pd.DataFrame(b); B.to_csv(OUT + "A4b_AAH_per_patient.tsv", sep="\t", index=False)
report["A4b"] = {"definition": "严格巨噬细胞（A1a 筛选）每病灶 SPP1 伪批量值（§6.1）；AAH 患者的 AAH、同患者浸润与正常病灶并列；仅描述，不做推断",
                 "n_patients_with_AAH": len(aah_pat), "patients": aah_pat,
                 "n_AAH_lesions_ge50_mac": int(((L.stage == "AAH") & L.qualifies).sum())}

# ======================= A4c：只差权重（§6.2） =======================
Q = L[L.qualifies & L.group.isin(["Precursor", "LUAD"])]
pats = sorted(set(Q[Q.group == "Precursor"].patient) & set(Q[Q.group == "LUAD"].patient), key=lambda x: int(x[1:]))
rows = []
for p in pats:
    rec = dict(patient=p)
    for g in ["Precursor", "LUAD"]:
        q = Q[(Q.patient == p) & (Q.group == g)]
        x, n = q.x.values, q.n_mac.values.astype(float)
        rec[f"k_{g}"] = int(len(q))
        rec[f"{g}_equal"] = float(x.sum() / len(x))                     # 版本一：病灶等权
        rec[f"{g}_cellweighted"] = float((n * x).sum() / n.sum())        # 版本二：细胞数加权（n_j = 合格巨噬细胞数）
        rec[f"{g}_lesions"] = ",".join(q.lesion)
    rows.append(rec)
W = pd.DataFrame(rows)
# 断言：同组只有一个合格病灶时两版必须完全相等（两版只差权重的运行时检查）
for g in ["Precursor", "LUAD"]:
    one = W[f"k_{g}"] == 1
    assert np.allclose(W.loc[one, f"{g}_equal"], W.loc[one, f"{g}_cellweighted"], rtol=0, atol=1e-12), g
res = {}
for ver in ["equal", "cellweighted"]:
    d = (W[f"LUAD_{ver}"] - W[f"Precursor_{ver}"]).values
    W[f"diff_{ver}"] = d
    t = stats.wilcoxon(W[f"LUAD_{ver}"], W[f"Precursor_{ver}"], alternative="greater")
    lo, hi = boot_ci(d, n=10000, seed=SEED)
    res[ver] = {"n_patients": int(len(d)), "n_positive": int((d > 0).sum()),
                "mean_diff": float(d.mean()), "median_diff": float(np.median(d)),
                "boot95CI_mean_diff": [lo, hi], "wilcoxon_one_sided_p": float(t.pvalue), "wilcoxon_stat": float(t.statistic)}
W.to_csv(OUT + "A4c_patient_values.tsv", sep="\t", index=False)
multi = W[(W.k_Precursor > 1) | (W.k_LUAD > 1)].patient.tolist()
report["A4c"] = {
    "qualifying_lesions": "严格巨噬细胞（Macrophage 且 EPCAM<2，与 A1a 同一筛选）≥50 个的病灶；两版共用同一病灶集合与同一患者集合",
    "x_j": "每合格病灶唯一的 x_j = log(1 + CP10k_j)，CP10k_j 由该病灶合格巨噬细胞原始计数聚合后以全基因总计数归一化（§6.1）",
    "version_equal": "患者内同组 Σx_j / k", "version_cellweighted": "患者内同组 Σ n_j·x_j / Σ n_j，n_j 为合格巨噬细胞数",
    "test": "配对 Wilcoxon 符号秩，单侧（LUAD > Precursor），scipy 默认零差处理；自助法 10000 次、种子 20261019，CI 针对患者级差值均值",
    "estimands": {"equal": "每患者病灶平均的差异", "cellweighted": "以巨噬细胞数加权的每患者差异"},
    "patients_where_weighting_can_differ": multi,
    "note": "两版只差权重；只有在同组含多于一个合格病灶的患者中两版数值才可能不同。A1a 原始计数负二项模型在 A1 中单独报告，不与本项比较。",
    "results": res}
log("A4c", {v: (r["n_patients"], round(r["mean_diff"], 3), r["wilcoxon_one_sided_p"]) for v, r in res.items()})
dump(report, OUT + "A4_results.json")
log("A4 完成")
