"""A1（方案 v1.3 §2-A1、§6.1、§6.4）严格巨噬细胞注释内的配对比较、分解与亚群敏感性。环境 ggn_target。
A1a  Macrophage 且 EPCAM<2；患者×病变组（前体 = AAH+AIS+MIA 合并，按细胞聚合计数）；每单位 ≥50 细胞；
     pydeseq2 ~ patient + condition；LUAD vs Precursor；SPP1 的 log2FC、SE、FDR；按方案三分支给出措辞类别。
A1b  同一组单位内分解：(i) SPP1+ 细胞比例；(ii) SPP1+ 细胞内 SPP1 强度（逐细胞 log1p(CP10k) 的均值，分母为该细胞全基因
     总计数）；单位内 SPP1+ <10 个时 (ii) 不可评估。患者级配对 Wilcoxon（方案未写方向 → 双侧），Holm 固定 m=2。
A1c  亚群敏感性（§6.4）：原始计数完整保留；归一化因子按全基因；待验证基因只从高变基因中剔除；
     PCA 50 → harmony（批次=样本）→ Leiden 0.5，种子 20261019；驻留样/募集样评分 score_genes 默认参数、random_state=0，
     在全部严格巨噬细胞上计算一次，群均值层面 z 标准化，(驻留−募集) ≥0.5 / ≤−0.5 / 其间 = Mixed/Unresolved；
     同表型多群先合并再检验；可配对患者 <8 记为不可评估；Holm m=2（驻留样 + 募集样）。
     亚群内检验按执行补充记录（docs/执行补充_A1c_A1d_2026-10-19.md，徐老师审定）：重复 A1a 模型，双侧；
     目标检验 = 两表型 SPP1 原始双侧 Wald P 的 Holm（固定 m=2）；全基因 BH padj 分开报告，不参与目标检验。
A1d（Visium 标记集敏感性）另见 A1d_visium_markers.py。"""
import os, sys, json, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump, holm_fixed, boot_ci, SEED
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
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
OUT = outdir("review_A1")
MIN_CELLS = 50; MIN_PAIRS = 8; MIN_POS = 10
PRE = ("AAH", "AIS", "MIA"); INV = "LUAD"
RESIDENT = ["MARCO", "FABP4", "MCEMP1", "PPARG"]
RECRUITED = ["VCAN", "FCN1", "S100A8", "S100A9"]

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

def paired_units(obs, mask, label):
    """患者×病变组单位（前体合并），两侧均 ≥50 细胞的患者才纳入。返回 {unit: 行号数组}, meta"""
    groups, meta, excl = {}, [], []
    for p in sorted(obs.patient.unique(), key=lambda x: int(x[1:])):
        idx = {g: np.where(mask & (obs.patient.values == p) & (obs.grp.values == g))[0] for g in ["Precursor", "LUAD"]}
        if min(len(idx["Precursor"]), len(idx["LUAD"])) < MIN_CELLS:
            excl.append(p); continue
        for g in ["Precursor", "LUAD"]:
            groups[f"{p}_{g}"] = idx[g]; meta.append(dict(unit=f"{p}_{g}", patient=p, condition=g, n_cells=len(idx[g]), subset=label))
    return groups, (pd.DataFrame(meta).set_index("unit") if meta else pd.DataFrame()), excl

def run_deseq(C, genes, groups, meta, save=None):
    counts = pd.DataFrame({k: np.asarray(C[i].sum(axis=0)).ravel() for k, i in groups.items()}, index=genes).T.astype(int)
    m = meta.copy(); m["condition"] = pd.Categorical(m["condition"], categories=["Precursor", "LUAD"])
    r = deseq(counts, m[["patient", "condition"]], "~ patient + condition", ["condition", "LUAD", "Precursor"])
    if save: r.to_csv(save, sep="\t")                    # 全基因结果单独保存
    s = r.loc["SPP1"]
    return dict(SPP1_log2FC=float(s.log2FoldChange), SPP1_lfcSE=float(s.lfcSE),
                SPP1_p_two_sided=float(s.pvalue), SPP1_padj=float(s.padj))

# ======================= 载入 =======================
a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad"); log("载入", a.shape)
assert "counts" in a.layers, "缺少 layers['counts']"
Call = a.layers["counts"].tocsr(); genes = a.var_names.copy(); gi = {g: i for i, g in enumerate(genes)}
ep = np.asarray(Call[:, gi["EPCAM"]].todense()).ravel()
mac = (a.obs["final_label"].astype(str).values == "Macrophage") & (ep < 2)            # A1a 严格注释
b = a[mac].copy(); del a, Call
b.X = b.layers["counts"].copy()
ob = b.obs[["patient", "sample", "stage"]].astype(str).copy()
ob["grp"] = np.where(ob.stage == INV, "LUAD", np.where(ob.stage.isin(PRE), "Precursor", "Other"))
C = b.layers["counts"].tocsr()
spp1 = np.asarray(C[:, gi["SPP1"]].todense()).ravel()
lib = np.asarray(C.sum(axis=1)).ravel()
log("严格巨噬细胞", b.n_obs)
report = {"plan_sha256": sha, "n_strict_macrophages": int(b.n_obs)}
allmask = np.ones(b.n_obs, bool)

# ======================= A1a =======================
groups, meta, excl = paired_units(ob, allmask, "all_strict_macrophages")
a1a = {"definition": "Macrophage 且 EPCAM<2；前体 = AAH+AIS+MIA 按细胞合并；每单位 ≥50；pydeseq2 ~ patient + condition；LUAD vs Precursor",
       "n_pairs": len(meta) // 2, "patients_excluded_lt50": excl}
a1a.update(run_deseq(C, genes, groups, meta, save=OUT + "A1a_DE_allgenes.tsv"))
lfc, fdr = a1a["SPP1_log2FC"], a1a["SPP1_padj"]
if lfc > 0 and fdr < 0.05:
    a1a["decision_branch"] = "方向一致且 FDR<0.05：可写'在所采用的巨噬细胞注释内，SPP1 表达差异仍然存在'；不得写成已排除组成变化或笼统'状态改变'"
elif lfc > 0:
    a1a["decision_branch"] = "方向一致但不显著：写'方向一致但在严格注释内未达显著'"
else:
    a1a["decision_branch"] = "方向反转：主张收窄为'髓系区室内 SPP1 升高'，删去全部'巨噬细胞状态改变'表述（含摘要）"
meta.to_csv(OUT + "A1a_units.tsv", sep="\t")
report["A1a"] = a1a; log("A1a", {k: a1a[k] for k in ["n_pairs", "SPP1_log2FC", "SPP1_padj"]})

# ======================= A1b =======================
rows = []
for u, idx in groups.items():
    pos = idx[spp1[idx] > 0]
    inten = float(np.mean(np.log1p(spp1[pos] / lib[pos] * 1e4))) if len(pos) >= MIN_POS else np.nan
    rows.append(dict(unit=u, patient=meta.loc[u, "patient"], condition=meta.loc[u, "condition"], n_cells=len(idx),
                     n_SPP1pos=int(len(pos)), frac_SPP1pos=float(len(pos) / len(idx)), intensity_SPP1pos=inten))
U = pd.DataFrame(rows); U.to_csv(OUT + "A1b_units.tsv", sep="\t", index=False)
a1b, pv = {}, {}
for k, col in [("i_fraction", "frac_SPP1pos"), ("ii_intensity", "intensity_SPP1pos")]:
    w = U.pivot(index="patient", columns="condition", values=col)
    n_unevaluable_units = int(U[col].isna().sum())
    w = w.dropna()
    d = (w["LUAD"] - w["Precursor"]).values
    if len(d) == 0 or np.all(d == 0):
        a1b[k] = {"evaluable": False, "n_patients": int(len(d)), "n_units_unevaluable": n_unevaluable_units}; pv[k] = None; continue
    t = stats.wilcoxon(w["LUAD"], w["Precursor"], alternative="two-sided")
    lo, hi = boot_ci(d, n=10000, seed=SEED)
    a1b[k] = {"evaluable": True, "n_patients": int(len(d)), "n_units_unevaluable": n_unevaluable_units,
              "n_LUAD_higher": int((d > 0).sum()), "mean_diff": float(d.mean()), "median_diff": float(np.median(d)),
              "boot95CI_mean_diff": [lo, hi], "wilcoxon_two_sided_p": float(t.pvalue)}
    pv[k] = float(t.pvalue)
for k, v in holm_fixed(pv, 2).items(): a1b[k]["holm_m2"] = v
a1b["note"] = "方案未写方向 → 双侧。(ii) 单位内 SPP1+ <10 记为不可评估；患者任一侧不可评估则不入 (ii) 的配对。"
report["A1b"] = a1b; log("A1b", {k: (v.get("n_patients"), v.get("holm_m2")) for k, v in a1b.items() if isinstance(v, dict)})

# ======================= A1c（§6.4） =======================
niche = pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv", sep="\t").gene.tolist()
gs = json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
HIF13 = ["CA9", "VEGFA", "SLC2A1", "PGK1", "LDHA", "BNIP3", "NDRG1", "ADM", "ANKRD37", "P4HA1", "ALDOA", "PDK1", "EGLN3"]
# HALLMARK_HYPOXIA 全部基因 = 106 实际所用（已去掉与 30 基因重叠者）∪ 30 基因中的重叠者；来源与发现队列一致（D1）
EXCL = sorted(set(["SPP1"]) | set(niche) | set(HIF13) | set(gs["HIF14_used"]) | set(gs["HALLMARK_used"]))
c = b.copy()
sc.pp.normalize_total(c, target_sum=1e4)           # 分母 = 全基因总计数（§6.4 第 2 条）
sc.pp.log1p(c)
c.layers["lognorm"] = c.X.copy()
sc.pp.highly_variable_genes(c, n_top_genes=2000)
hvg = c.var_names[c.var.highly_variable].tolist()
removed = [g for g in hvg if g in EXCL]
feat = [g for g in hvg if g not in EXCL]             # 只从聚类特征集合中排除（§6.4 第 3 条）
log("HVG", len(hvg), "剔除", len(removed), "保留", len(feat))
f = c[:, feat].copy()
sc.pp.scale(f, max_value=10)
sc.tl.pca(f, n_comps=50, random_state=SEED)
# 直接调用 harmonypy（scanpy 封装与 harmonypy ≥2 的返回方向不兼容）；按形状判定方向
import harmonypy as hm
ho = hm.run_harmony(f.obsm["X_pca"], f.obs[["sample"]].astype(str), "sample", random_state=SEED)
Z = np.asarray(ho.Z_corr)
Z = Z if Z.shape[0] == f.n_obs else Z.T
assert Z.shape == (f.n_obs, f.obsm["X_pca"].shape[1]), Z.shape
f.obsm["X_pca_harmony"] = Z
sc.pp.neighbors(f, use_rep="X_pca_harmony", random_state=SEED)
sc.tl.leiden(f, resolution=0.5, random_state=SEED, key_added="leiden_A1c")
cl = f.obs["leiden_A1c"].astype(str).values
del f
# 表型评分：全部严格巨噬细胞上一次（完整基因、对数归一化矩阵）
sc.tl.score_genes(c, [g for g in RESIDENT if g in c.var_names], score_name="resident", random_state=0)
sc.tl.score_genes(c, [g for g in RECRUITED if g in c.var_names], score_name="recruited", random_state=0)
CL = pd.DataFrame({"cluster": cl, "resident": c.obs["resident"].values, "recruited": c.obs["recruited"].values})
G = CL.groupby("cluster").agg(n_cells=("resident", "size"), resident_mean=("resident", "mean"), recruited_mean=("recruited", "mean"))
z = lambda s: (s - s.mean()) / s.std(ddof=0)
G["resident_z"] = z(G.resident_mean); G["recruited_z"] = z(G.recruited_mean)
G["diff_z"] = G.resident_z - G.recruited_z
G["phenotype"] = np.where(G.diff_z >= 0.5, "resident-like", np.where(G.diff_z <= -0.5, "monocyte-recruited-like", "Mixed/Unresolved"))
G["cell_frac"] = G.n_cells / G.n_cells.sum()
G.to_csv(OUT + "A1c_clusters.tsv", sep="\t")
del c
ph = G.phenotype.reindex(cl).values
a1c, pv = {"excluded_genes_from_HVG": removed, "n_excluded_genes_considered": len(EXCL), "n_HVG_used": len(feat),
           "test_rule": "执行补充记录 §1：亚群内重复 A1a 模型（pydeseq2 ~ patient + condition，单位 ≥50，完整患者对，<8 对不可评估），"
                        "双侧；目标检验 = 两表型 SPP1 原始双侧 Wald P 的 Holm（固定 m=2）。全基因 BH padj 另列，不参与目标检验。",
           "supplementary_record_sha256": _SUPP_SHA}, {}
for phen in ["resident-like", "monocyte-recruited-like"]:
    msk = (ph == phen)
    if msk.sum() == 0:
        a1c[phen] = {"evaluable": False, "reason": "无该表型的群", "n_cells": 0}; pv[phen] = None; continue
    g2, m2, ex2 = paired_units(ob, msk, phen)
    npairs = len(m2) // 2 if len(m2) else 0
    if npairs < MIN_PAIRS:
        a1c[phen] = {"evaluable": False, "reason": f"可配对患者 {npairs} < {MIN_PAIRS}", "n_cells": int(msk.sum()),
                     "clusters": G.index[G.phenotype == phen].tolist()}; pv[phen] = None; continue
    r = run_deseq(C, genes, g2, m2, save=OUT + f"A1c_DE_allgenes_{phen}.tsv")
    a1c[phen] = {"evaluable": True, "n_pairs": npairs, "n_cells": int(msk.sum()), "clusters": G.index[G.phenotype == phen].tolist(),
                 "SPP1_log2FC": r["SPP1_log2FC"], "SPP1_lfcSE": r["SPP1_lfcSE"], "SPP1_p_two_sided": r["SPP1_p_two_sided"],
                 "genomewide_BH_padj_SPP1_not_used_for_target_test": r["SPP1_padj"]}
    pv[phen] = r["SPP1_p_two_sided"]                       # 目标检验：原始双侧 P，不用全基因 padj
    m2.to_csv(OUT + f"A1c_units_{phen}.tsv", sep="\t")
for k, v in holm_fixed(pv, 2).items(): a1c[k]["target_test_holm_m2_on_raw_two_sided_p"] = v
mixed = G[G.phenotype == "Mixed/Unresolved"]
a1c["Mixed/Unresolved"] = {"n_clusters": int(len(mixed)), "n_cells": int(mixed.n_cells.sum()), "cell_frac": float(mixed.cell_frac.sum()),
                           "note": "不做分组推断，仅报告"}
a1c["note"] = "表型名称仅为表达表型（resident-like / monocyte-recruited-like），不表示已证明的细胞来源；A1c 为敏感性分析，不替代 A1a。"
report["A1c"] = a1c
log("A1c", {k: (v.get("n_pairs"), v.get("target_test_holm_m2_on_raw_two_sided_p")) for k, v in a1c.items() if isinstance(v, dict) and "evaluable" in v})
dump(report, OUT + "A1_results.json")
log("A1 完成")
