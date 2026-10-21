"""A6（方案 v1.3 §2-A6、§6.5）预后模型的完整报告 —— 第 1 步：Python 部分。环境 ggn_target。
复现原报告的全部队列内 Cox 模型（statsmodels PHReg，Breslow；与 96/98/102/105 逐项相同），与表 S6b 核对；
记录分析人群（分期合格数、因随访时间缺失/≤0 剔除数、无表达数据剔除数、年龄插补数、性别原始取值与编码）；
导出每个队列内模型的设计矩阵供 R 做 Schoenfeld 诊断与 Firth（A6_cox_diagnostics.R）；
合并分析（DerSimonian–Laird、Hartung–Knapp、τ²、I²、95% PI）含留一与分终点；A6d 共线性；A6e 基因覆盖；A6f 基因来源。
不修改 09_layer4_clinical.py：载入时对其源码副本插入记录语句后执行，逻辑不变。"""
import os, sys, json, types, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
warnings.filterwarnings("ignore")

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A6"); os.makedirs(OUT + "designs", exist_ok=True)
sys.path.insert(0, f"{P}/scripts")
from prolif_genes import TIROSH

# ---------------- 载入 09（源码副本插入记录语句，逻辑不变） ----------------
REC = {}
src_path = f"{P}/scripts/09_layer4_clinical.py"
src = open(src_path, encoding="utf-8").read()
INJ = [
 ('    d = d.dropna(subset=["time", "event"]); d = d[d["time"] > 0]\n',
  '    _REC[name] = {"n_stage_eligible": int(len(d)), "age_na_ids": list(d.index[d["age"].isna()]), "clin": c}\n'
  '    d = d.dropna(subset=["time", "event"]); d = d[d["time"] > 0]\n'
  '    _REC[name]["n_after_time_filter"] = int(len(d))\n'),
 ('    return d.loc[g.columns], g, plat, endpoint, raw_scale\n',
  '    _REC[name]["n_with_expression"] = int(len(g.columns))\n'
  '    return d.loc[g.columns], g, plat, endpoint, raw_scale\n'),
 ('    d = d.dropna(subset=["time"]); d = d[d["time"] > 0]; d["age"] = d["age"].fillna(d["age"].median())\n',
  '    _REC["TCGA"] = {"n_stage_eligible": int(len(d)), "age_na_ids": list(d.index[d["age"].isna()]), "clin": c}\n'
  '    d = d.dropna(subset=["time"]); d = d[d["time"] > 0]; d["age"] = d["age"].fillna(d["age"].median())\n'
  '    _REC["TCGA"]["n_after_time_filter"] = int(len(d))\n'),
 ('    return d.loc[cols], e, "RNA-seq", "OS", False\n',
  '    _REC["TCGA"]["n_with_expression"] = int(len(cols))\n'
  '    return d.loc[cols], e, "RNA-seq", "OS", False\n'),
]
for o, n in INJ:
    assert src.count(o) == 1, f"09 源码锚点不唯一或不存在：{o[:60]}"
    src = src.replace(o, n)
l4 = types.ModuleType("l4"); l4.__dict__["_REC"] = REC; l4.__file__ = src_path
exec(compile(src, src_path, "exec"), l4.__dict__)

NAMES = ["GSE31210", "GSE50081", "GSE30219", "GSE37745", "GSE68465", "GSE72094", "TCGA"]
cohorts = {n: (l4.tcga() if n == "TCGA" else l4.cohort_def(n)) for n in NAMES}
log("队列载入", {n: len(c[0]) for n, c in cohorts.items()})

# ---------------- 基因集（与原脚本逐字相同） ----------------
MYE = "CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()
CAF = ["FAP", "POSTN", "CTHRC1", "COL11A1", "COL10A1", "MMP11", "INHBA"]
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]; IMM = ["PTPRC"]
HYP14 = ["CA9", "VEGFA", "SLC2A1", "PGK1", "LDHA", "BNIP3", "NDRG1", "ADM", "ANKRD37", "P4HA1", "ENO1", "ALDOA", "PDK1", "EGLN3"]  # 102 实际所用（含 ENO1）
HIF13 = [g for g in HYP14 if g != "ENO1"]                                                                                         # 稿件所称 HIF13
top = pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv", sep="\t").gene.tolist(); assert len(top) == 30
rat = pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv", sep="\t", index_col=0)
best = rat[["epi"] + [c for c in rat.columns if c.startswith("med_")]].idxmax(axis=1)
top_mye = [g for g in top if not (g in best.index and best[g] != "med_Myeloid")]             # 同 98
try:
    import gseapy as gp
    HALL = gp.get_library("MSigDB_Hallmark_2020").get("Hypoxia"); assert HALL and len(HALL) > 100
except Exception as e:
    HALL = None; log("HALLMARK 取不到（+HALLMARK 模型记为无法复现）：", repr(e))
GS = {"mye": MYE, "caf": CAF, "epi": EPI, "imm": IMM, "hyp14": HYP14, "hif13": HIF13, "prolif": TIROSH, "hall": HALL}

def zr(m): return m.sub(m.mean(1), axis=0).div(m.std(1) + 1e-9, axis=0)
def score(g, genes):
    pres = [x for x in genes if x in g.index]; s = zr(g.loc[pres].astype(float)).mean(0)
    return (s - s.mean()) / (s.std() + 1e-9), len(pres)

MODELS = {   # 名称: (评分基因, 调整项, 稿件行名, 来源脚本)
 "main":               (top, ["mye"], "Main model", "96"),
 "plus_prolif":        (top, ["mye", "prolif"], "+ proliferation score", "96"),
 "without_SPP1":       ([x for x in top if x != "SPP1"], ["mye"], "Without SPP1 in score", "96"),
 "SPP1_alone":         (["SPP1"], ["mye"], "SPP1 alone", "96"),
 "unadjusted_myeloid": (top, [], "Not adjusted for myeloid score", "96"),
 "plus_CAF":           (top, ["mye", "caf"], "+ CAF score", "98"),
 "myeloid15":          (top_mye, ["mye"], "Myeloid-source genes only (15)", "98"),
 "plus_purity":        (top, ["mye", "epi", "imm"], "+ purity proxy", "102"),
 "plus_HIF14_reported_as_HIF13": (top, ["mye", "hyp14"], "+ hypoxia (HIF13)", "102"),
 "plus_HALLMARK":      (top, ["mye", "hall"], "+ hypoxia (HALLMARK)", "105"),
 "fully_adjusted":     (top, ["mye", "epi", "imm", "hyp14", "prolif", "caf"], "Fully adjusted", "102"),
 # 更正候选（方案 0.3）：稿件标为 HIF13，原代码用了含 ENO1 的 14 基因
 "CORR_plus_HIF13":       (top, ["mye", "hif13"], "(correction candidate) + hypoxia, HIF13 without ENO1", "A6"),
 "CORR_fully_adjusted_HIF13": (top, ["mye", "epi", "imm", "hif13", "prolif", "caf"], "(correction candidate) fully adjusted, HIF13", "A6"),
}
EXPECTED = {  # 表 S6b 已报告值（两位小数），用于复现核对
 "main": (1.59, 1.11, 2.29, 0.71), "plus_prolif": (1.44, 1.02, 2.05, 0.65), "without_SPP1": (1.59, 1.09, 2.30, 0.72),
 "SPP1_alone": (1.31, 1.08, 1.59, 0.47), "unadjusted_myeloid": (1.20, 0.97, 1.49, 0.67), "plus_CAF": (1.64, 1.06, 2.52, 0.71),
 "myeloid15": (1.65, 1.17, 2.33, 0.66), "plus_purity": (1.63, 1.11, 2.41, 0.71), "plus_HIF14_reported_as_HIF13": (1.26, 0.96, 1.65, 0.30),
 "plus_HALLMARK": (1.43, 1.01, 2.03, 0.57), "fully_adjusted": (1.53, 0.95, 2.44, 0.66),
 "LOO_GSE31210": (1.36, 1.04, 1.79, 0.37), "LOO_GSE50081": (1.58, 1.05, 2.37, 0.75), "LOO_GSE30219": (1.53, 1.03, 2.27, 0.74),
 "LOO_GSE37745": (1.68, 1.15, 2.45, 0.74), "LOO_GSE68465": (1.80, 1.28, 2.53, 0.52), "LOO_GSE72094": (1.55, 1.00, 2.39, 0.74),
 "LOO_TCGA": (1.66, 1.06, 2.59, 0.75), "endpoint_RFS_DFS": (1.61, 0.91, 2.86, 0.79), "endpoint_OS": (1.53, 1.10, 2.12, 0.18),
}

# ---------------- 队列内 Cox（与原脚本相同的设计） ----------------
def design(d, s, extra):
    X0 = d[["age", "male", "stage"]].copy(); X0["age"] = (X0.age - X0.age.mean()) / (X0.age.std() + 1e-9)
    dropped = [c for c in X0.columns if not X0[c].std() > 0]; X0 = X0.loc[:, X0.std() > 0]
    for k, v in extra.items(): X0[k] = v.reindex(d.index).values
    X = pd.concat([s.reindex(d.index).rename("score"), X0], axis=1)
    return X, dropped

rows, designs = [], []
for mname, (gs, adj, label, origin) in MODELS.items():
    if "hall" in adj and HALL is None:
        rows.append(dict(model=mname, label=label, origin=origin, cohort="ALL", fit_status="not reproducible: HALLMARK library unavailable")); continue
    for n, c in cohorts.items():
        d, g, plat, ep = c[0], c[1], c[2], c[3]
        s, np_ = score(g, gs); extra, cov_n = {}, {}
        for a_ in adj: extra[a_], cov_n[a_] = score(g, GS[a_])
        X, dropped = design(d, s, extra)
        rec = dict(model=mname, label=label, origin=origin, cohort=n, endpoint=ep, n=int(len(d)), events=int(d.event.sum()),
                   n_covariates_incl_score=int(X.shape[1]), df=int(X.shape[1]), covariates=";".join(X.columns),
                   constant_covariates_dropped=";".join(dropped), score_genes_present=int(np_), score_genes_total=int(len(gs)),
                   covariate_genes_present=json.dumps(cov_n))
        try:
            r = PHReg(d.time.values, X.values, status=d.event.values, ties="breslow").fit(disp=False)
            conv = r.mle_retvals.get("converged", None) if hasattr(r, "mle_retvals") and r.mle_retvals else None
            b, e = float(r.params[0]), float(r.bse[0])
            rec.update(beta=b, se=e, HR=float(np.exp(b)), lo=float(np.exp(b - 1.96 * e)), hi=float(np.exp(b + 1.96 * e)),
                       p_one_sided=float(norm.sf(b / e)), p_two_sided=float(2 * norm.sf(abs(b / e))),
                       converged=conv, fit_status="ok" if np.isfinite(b) and np.isfinite(e) else "non-finite")
        except Exception as ex:
            rec.update(fit_status=f"failed: {ex!r}")
        rows.append(rec)
        fn = f"{mname}__{n}.tsv"
        pd.concat([d[["time", "event"]], X], axis=1).to_csv(OUT + "designs/" + fn, sep="\t", index_label="sample")
        designs.append(dict(model=mname, cohort=n, file=fn, events=int(d.event.sum())))
M = pd.DataFrame(rows); M.to_csv(OUT + "A6_cohort_models_python.tsv", sep="\t", index=False)
pd.DataFrame(designs).to_csv(OUT + "designs/manifest.tsv", sep="\t", index=False)

# ---------------- 合并分析（§6.5(b)）----------------
def meta(lh, se, cohorts_used):
    lh, se = np.array(lh), np.array(se); k = len(lh); w = 1 / se ** 2; fe = (w * lh).sum() / w.sum(); Q = (w * (lh - fe) ** 2).sum()
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); w2 = 1 / (se ** 2 + tau2); mu = (w2 * lh).sum() / w2.sum(); s = np.sqrt(1 / w2.sum())
    pi = stats.t.ppf(0.975, k - 2) * np.sqrt(tau2 + s ** 2) if k > 2 else np.nan
    qhk = (w2 * (lh - mu) ** 2).sum() / (k - 1); shk = np.sqrt(qhk / w2.sum()); thk = stats.t.ppf(0.975, k - 1)
    return dict(k=k, cohorts=";".join(cohorts_used), method="DerSimonian-Laird random effects",
                HR=float(np.exp(mu)), lo=float(np.exp(mu - 1.96 * s)), hi=float(np.exp(mu + 1.96 * s)),
                p_one_sided=float(norm.sf(mu / s)), p_two_sided=float(2 * norm.sf(abs(mu / s))),
                tau2=float(tau2), I2=float(max(0, (Q - (k - 1)) / Q)) if Q > 0 else 0.0,
                PI95_lo=float(np.exp(mu - pi)) if k > 2 else np.nan, PI95_hi=float(np.exp(mu + pi)) if k > 2 else np.nan,
                HK_lo=float(np.exp(mu - thk * shk)), HK_hi=float(np.exp(mu + thk * shk)),
                HK_p_two_sided=float(2 * stats.t.sf(abs(mu / shk), k - 1)))
pooled = []
ok = M[M.get("fit_status", pd.Series(dtype=str)) == "ok"]
for mname in MODELS:
    m = ok[ok.model == mname]
    if len(m) == len(NAMES): pooled.append(dict(model=mname, **meta(m.beta, m.se, list(m.cohort))))
    else: pooled.append(dict(model=mname, k=len(m), method="not pooled: not all cohorts fitted"))
mm = ok[ok.model == "main"].set_index("cohort")
for n in NAMES:
    sub = mm.drop(index=n); pooled.append(dict(model=f"LOO_{n}", **meta(sub.beta, sub.se, list(sub.index))))
for epn, eps in [("endpoint_RFS_DFS", ["RFS", "DFS"]), ("endpoint_OS", ["OS"])]:
    sub = mm[mm.endpoint.isin(eps)]
    rec = dict(model=epn, **meta(sub.beta, sub.se, list(sub.index))); rec.update(n_patients=int(sub.n.sum()), n_events=int(sub.events.sum()))
    pooled.append(rec)
Pm = pd.DataFrame(pooled)
chk = []
for mname, (h, l, u, i2) in EXPECTED.items():
    r = Pm.set_index("model").loc[mname] if mname in set(Pm.model) else None
    if r is None or "HR" not in r or pd.isna(r.get("HR", np.nan)):
        chk.append(dict(model=mname, reproduced=False, note="not computed")); continue
    got = (round(r.HR, 2), round(r.lo, 2), round(r.hi, 2), round(r.I2, 2))
    chk.append(dict(model=mname, reported=str((h, l, u, i2)), recomputed=str(got), reproduced=bool(got == (h, l, u, i2))))
C = pd.DataFrame(chk); Pm = Pm.merge(C[["model", "reproduced"]], on="model", how="left")
Pm.to_csv(OUT + "A6_pooled.tsv", sep="\t", index=False); C.to_csv(OUT + "A6_reproduction_check.tsv", sep="\t", index=False)
log("复现核对", int(C.reproduced.sum()), "/", len(C))

# ---------------- 分析人群与缺失值处理 ----------------
SEXCOL = ["gender", "Sex", "sex", "gender.demographic"]
pop = []
for n, c in cohorts.items():
    d = c[0]; R = REC[n]; clin = R["clin"]
    sc_ = next((k for k in SEXCOL if k in clin.columns), None)
    raw = clin.loc[d.index, sc_].astype(str).str.strip() if sc_ else pd.Series(dtype=str)
    vc = raw.value_counts(dropna=False).to_dict()
    recognised = raw.str.lower().isin(["male", "m", "female", "f"])
    pop.append(dict(cohort=n, platform=c[2], endpoint=c[3], n_stage_eligible=R["n_stage_eligible"],
                    excluded_missing_or_nonpositive_time=R["n_stage_eligible"] - R["n_after_time_filter"],
                    excluded_no_expression=R["n_after_time_filter"] - R["n_with_expression"],
                    n_analysed=int(len(d)), events=int(d.event.sum()),
                    age_missing_median_imputed=int(len(set(R["age_na_ids"]) & set(d.index))),
                    sex_column=sc_, sex_raw_values=json.dumps({str(k): int(v) for k, v in vc.items()}, ensure_ascii=False),
                    sex_unrecognised_coded_as_0=int((~recognised).sum()) if sc_ else None,
                    n_male_coded=int(d.male.sum()), n_stage_IB_or_T2=int(d.stage.sum())))
POP = pd.DataFrame(pop); POP.to_csv(OUT + "A6_analysis_population.tsv", sep="\t", index=False)

# ---------------- A6d 共线性（主模型人群；调整前后模型人群相同） ----------------
col = []
for n, c in cohorts.items():
    d, g = c[0], c[1]
    s = score(g, top)[0]; m_ = score(g, MYE)[0]
    X, _ = design(d, s, {"mye": m_})
    def vif(col_):
        y = X[col_].values; Z = np.column_stack([np.ones(len(X)), X.drop(columns=col_).values])
        r2 = 1 - ((y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        return float(1 / (1 - r2))
    col.append(dict(cohort=n, n=int(len(d)), spearman_score_vs_myeloid=float(stats.spearmanr(s, m_)[0]),
                    pearson_score_vs_myeloid=float(stats.pearsonr(s, m_)[0]), VIF_score=vif("score"), VIF_myeloid=vif("mye"),
                    same_population_adjusted_vs_unadjusted=True))
COL = pd.DataFrame(col); COL.to_csv(OUT + "A6_collinearity.tsv", sep="\t", index=False)

# ---------------- A6e 基因覆盖 ----------------
cov = []
for n, c in cohorts.items():
    g = c[1]; pres = [x for x in top if x in g.index]
    cov.append(dict(cohort=n, platform=c[2], genes_measured=len(pres), genes_missing=";".join(x for x in top if x not in g.index),
                    mapping=("Ensembl to symbol via GENCODE v44; duplicate symbols resolved by highest mean expression" if n == "TCGA" else
                             "probe to symbol via platform annotation; probes mapping to several genes (///) excluded; per gene, probe with highest mean expression"),
                    log2_transformed_from_MAS5=bool(c[4])))
pd.DataFrame(cov).to_csv(OUT + "A6_gene_coverage.tsv", sep="\t", index=False)

# ---------------- A6f 基因来源核查 ----------------
de = pd.read_csv(f"{P}/results/30_spatial/clinical/MacSPP1_vs_Macother_DE.tsv.gz", sep="\t", index_col=0)
IG = [x for x in top if x.startswith(("IGH", "IGK", "IGL", "JCHAIN"))]; COLG = [x for x in top if x.startswith("COL")]
src_rows = [dict(gene=x, family=("immunoglobulin" if x in IG else "collagen"),
                 detection_fraction_in_Mac_SPP1=float(de.loc[x, "det_SPP1grp"]) if x in de.index else np.nan,
                 myeloid_to_epithelial_median_ratio=float(de.loc[x, "src_ratio"]) if x in de.index else np.nan) for x in IG + COLG]
pd.DataFrame(src_rows).to_csv(OUT + "A6_gene_source.tsv", sep="\t", index=False)

dump({"plan_sha256": sha, "n_models_cohort_level": len(MODELS), "events_lt20_cohorts": [n for n, c in cohorts.items() if c[0].event.sum() < 20],
      "HALLMARK_available": HALL is not None, "HALLMARK_n_genes": len(HALL) if HALL else None,
      "reproduction": C.to_dict("records"),
      "hif_label_issue": "表 S6b '+ hypoxia (HIF13)' 与 'Fully adjusted' 由 102 用含 ENO1 的 14 基因计算；按原样复现并另算 HIF13（去 ENO1）更正候选，由研究者按 0.3 决定",
      "A6f_genes": {"immunoglobulin": IG, "collagen": COLG},
      "A6f_note": "检出比例只能提示来源疑问，既不能证明也不能排除其他细胞来源的污染（方案 A6f）"},
     OUT + "A6_prepare.json")
log("A6 第 1 步完成；下一步：Rscript scripts/review_A/A6_cox_diagnostics.R", OUT)
