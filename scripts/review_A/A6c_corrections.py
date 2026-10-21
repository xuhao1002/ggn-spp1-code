"""A6c（偏离记录 D4；方案 0.3 更正）批量队列中缺氧基因集按定义去除与评分重叠的基因后重算；GSE37745 剔除者比较。环境 ggn_target。
1) 重算：Cox "+HIF13"、"+HALLMARK"、"Fully adjusted"（队列内 + DL/HK/τ²/I²/PI）；评分与缺氧评分的逐队列 Spearman ρ 及中位数。
   每项同时按原代码复现原值（HIF 14 基因含 ENO1；HALLMARK 全部 200 基因），作对照与复现核对。
2) GSE37745：分期合格者中，因复发随访时间缺失/≤0 被剔除者 vs 保留者的年龄、性别、分期、总生存信息可得性（描述性；
   Mann–Whitney / Fisher 双侧 P 仅供描述）。
队列载入、评分函数、Cox 设计与原脚本 102/105 逐项相同。HALLMARK 来源为 Enrichr MSigDB_Hallmark_2020（D1）。"""
import os, sys, json, importlib.util, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
warnings.filterwarnings("ignore")

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A6c")
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("l4", f"{P}/scripts/09_layer4_clinical.py"); l4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
from prolif_genes import TIROSH
import gseapy as gp
HALL200 = gp.get_library("MSigDB_Hallmark_2020")["Hypoxia"]; assert len(HALL200) == 200, len(HALL200)

MYE = "CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()
CAF = ["FAP", "POSTN", "CTHRC1", "COL11A1", "COL10A1", "MMP11", "INHBA"]; EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]; IMM = ["PTPRC"]
HIF14 = ["CA9", "VEGFA", "SLC2A1", "PGK1", "LDHA", "BNIP3", "NDRG1", "ADM", "ANKRD37", "P4HA1", "ENO1", "ALDOA", "PDK1", "EGLN3"]
top = pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv", sep="\t").gene.tolist(); assert len(top) == 30
HIF13 = [g for g in HIF14 if g not in top]; HALL197 = [g for g in HALL200 if g not in top]
assert len(HIF13) == 13 and "ENO1" not in HIF13 and len(HALL197) == 197, (len(HIF13), len(HALL197))
gs_used = json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
assert set(HIF13) == set(gs_used["HIF14_used"]) and set(HALL197) == set(gs_used["HALLMARK_used"]), "与 106 的定义不一致"
log("HIF13 与 HALLMARK197 与 106 输出一致；重叠基因", sorted(set(top) & (set(HIF14) | set(HALL200))))

NAMES = ["GSE31210", "GSE50081", "GSE30219", "GSE37745", "GSE68465", "GSE72094", "TCGA"]
cohorts = {n: (l4.tcga() if n == "TCGA" else l4.cohort_def(n)) for n in NAMES}
def zr(m): return m.sub(m.mean(1), axis=0).div(m.std(1) + 1e-9, axis=0)
def score(g, genes):
    pres = [x for x in genes if x in g.index]; s = zr(g.loc[pres].astype(float)).mean(0); return (s - s.mean()) / (s.std() + 1e-9), len(pres)
def cox(d, x, extra):
    X0 = d[["age", "male", "stage"]].copy(); X0["age"] = (X0.age - X0.age.mean()) / (X0.age.std() + 1e-9); X0 = X0.loc[:, X0.std() > 0]
    for k, v in extra.items(): X0[k] = v.reindex(d.index).values
    r = PHReg(d.time.values, pd.concat([x.reindex(d.index).rename("s"), X0], axis=1).values, status=d.event.values, ties="breslow").fit(disp=False)
    return float(r.params[0]), float(r.bse[0])
def meta(lh, se):
    lh, se = np.array(lh), np.array(se); k = len(lh); w = 1 / se ** 2; fe = (w * lh).sum() / w.sum(); Q = (w * (lh - fe) ** 2).sum()
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); w2 = 1 / (se ** 2 + tau2); mu = (w2 * lh).sum() / w2.sum(); s = np.sqrt(1 / w2.sum())
    pi = stats.t.ppf(0.975, k - 2) * np.sqrt(tau2 + s ** 2)
    qhk = (w2 * (lh - mu) ** 2).sum() / (k - 1); shk = np.sqrt(qhk / w2.sum()); thk = stats.t.ppf(0.975, k - 1)
    return dict(HR=float(np.exp(mu)), lo=float(np.exp(mu - 1.96 * s)), hi=float(np.exp(mu + 1.96 * s)), p_one_sided=float(norm.sf(mu / s)),
                p_two_sided=float(2 * norm.sf(abs(mu / s))), tau2=float(tau2), I2=float(max(0, (Q - (k - 1)) / Q)) if Q > 0 else 0.0,
                PI95=[float(np.exp(mu - pi)), float(np.exp(mu + pi))], HK95=[float(np.exp(mu - thk * shk)), float(np.exp(mu + thk * shk))], k=k)
COV = {"mye": MYE, "epi": EPI, "imm": IMM, "prolif": TIROSH, "caf": CAF, "hif14": HIF14, "hif13": HIF13, "hall200": HALL200, "hall197": HALL197}
MODELS = {"plus_HIF_original_14": ["mye", "hif14"], "plus_HIF13_corrected": ["mye", "hif13"],
          "plus_HALLMARK_original_200": ["mye", "hall200"], "plus_HALLMARK197_corrected": ["mye", "hall197"],
          "fully_adjusted_original_HIF14": ["mye", "epi", "imm", "hif14", "prolif", "caf"],
          "fully_adjusted_corrected_HIF13": ["mye", "epi", "imm", "hif13", "prolif", "caf"]}
REPORTED = {"plus_HIF_original_14": (1.26, 0.96, 1.65), "plus_HALLMARK_original_200": (1.43, 1.01, 2.03), "fully_adjusted_original_HIF14": (1.53, 0.95, 2.44)}
res, per = {"plan_sha256": sha, "deviation": "D4"}, []
for m, adj in MODELS.items():
    lh, se = [], []
    for n, c in cohorts.items():
        d, g = c[0], c[1]; s, _ = score(g, top); ex, cn = {}, {}
        for a in adj: ex[a], cn[a] = score(g, COV[a])
        b, e = cox(d, s, ex); lh.append(b); se.append(e)
        per.append(dict(model=m, cohort=n, n=len(d), events=int(d.event.sum()), HR=float(np.exp(b)), lo=float(np.exp(b - 1.96 * e)),
                        hi=float(np.exp(b + 1.96 * e)), beta=b, se=e, covariate_genes_present=json.dumps(cn)))
    res[m] = meta(lh, se)
    if m in REPORTED:
        got = (round(res[m]["HR"], 2), round(res[m]["lo"], 2), round(res[m]["hi"], 2))
        res[m]["reported"] = REPORTED[m]; res[m]["reproduces_reported"] = bool(got == REPORTED[m])
    log(m, round(res[m]["HR"], 3), round(res[m]["lo"], 3), round(res[m]["hi"], 3))
pd.DataFrame(per).to_csv(OUT + "A6c_cox_per_cohort.tsv", sep="\t", index=False)

rho = {}
for lab, gs in [("HIF_original_14", HIF14), ("HIF13_corrected", HIF13), ("HALLMARK_original_200", HALL200), ("HALLMARK197_corrected", HALL197)]:
    rho[lab] = {n: float(stats.spearmanr(score(c[1], top)[0], score(c[1], gs)[0])[0]) for n, c in cohorts.items()}
res["rho_per_cohort"] = rho; res["rho_median"] = {k: float(np.median(list(v.values()))) for k, v in rho.items()}
res["rho_reported_median"] = {"HIF13 (as reported)": 0.35, "HALLMARK (as reported)": 0.58}
log("ρ 中位数", {k: round(v, 3) for k, v in res["rho_median"].items()})

# ---- GSE37745：被剔除 vs 保留 ----
clin, _ = l4.read_series(f"{P}/raw/stageI_cohorts/GSE37745_series_matrix.txt.gz")
c = clin[(clin["histology"] == "adeno") & clin["tumor stage"].isin(["1a", "1b"])].copy()
t = pd.to_numeric(c["days to recurrence / to last visit"], errors="coerce")
c["retained_candidate"] = t.notna() & (t > 0)
kept = set(cohorts["GSE37745"][0].index)
c["retained"] = c.index.isin(kept)
assert len(c) == 70 and c.retained.sum() == 34, (len(c), int(c.retained.sum()))
excl_expr = int((c.retained_candidate & ~c.retained).sum())
age = pd.to_numeric(c["age"], errors="coerce")
def summ(mask):
    s = c[mask]; a = age[mask]
    return dict(n=int(mask.sum()), age_median=float(a.median()), age_IQR=[float(a.quantile(.25)), float(a.quantile(.75))], age_missing=int(a.isna().sum()),
                male=int((s["gender"] == "male").sum()), female=int((s["gender"] == "female").sum()),
                stage_1a=int((s["tumor stage"] == "1a").sum()), stage_1b=int((s["tumor stage"] == "1b").sum()),
                recurrence_values=s["recurrence"].astype(str).value_counts().to_dict(),
                os_time_available=int(pd.to_numeric(s["days to determined death status"], errors="coerce").notna().sum()),
                dead_values=s["dead"].astype(str).value_counts().to_dict(),
                adjuvant_values=s["adjuvant treatment"].astype(str).value_counts().to_dict() if "adjuvant treatment" in s else None)
K, X = c.retained.values, (~c.retained).values
cmp = {"retained": summ(K), "excluded": summ(X), "excluded_with_time_but_no_expression": excl_expr}
cmp["age_mannwhitney_p_two_sided"] = float(stats.mannwhitneyu(age[K].dropna(), age[X].dropna(), alternative="two-sided").pvalue)
cmp["sex_fisher_p_two_sided"] = float(stats.fisher_exact([[cmp["retained"]["male"], cmp["retained"]["female"]], [cmp["excluded"]["male"], cmp["excluded"]["female"]]])[1])
cmp["stage_fisher_p_two_sided"] = float(stats.fisher_exact([[cmp["retained"]["stage_1a"], cmp["retained"]["stage_1b"]], [cmp["excluded"]["stage_1a"], cmp["excluded"]["stage_1b"]]])[1])
cmp["note"] = "描述性比较；P 值仅供描述，不作推断。剔除原因：复发随访时间（days to recurrence / to last visit）缺失或 ≤0。"
res["GSE37745_excluded_vs_retained"] = cmp
c[["age", "gender", "tumor stage", "recurrence", "days to recurrence / to last visit", "dead", "days to determined death status", "retained"]].to_csv(OUT + "A6c_GSE37745_eligible.tsv", sep="\t")
log("GSE37745", {k: cmp[k] for k in ["age_mannwhitney_p_two_sided", "sex_fisher_p_two_sided", "stage_fisher_p_two_sided"]})
dump(res, OUT + "A6c_results.json")
log("A6c 完成")
