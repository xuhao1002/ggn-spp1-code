"""偏离记录 D5（post hoc，2026-10-21 13:52 写于运行前）：GSE37745 总生存（OS）分析。环境 ggn_target。
三个拟合：(a) 全部 OS 可用者；(b) RFS 分析所含 34 例（终点换为 OS）；(c) 被排除的 36 例。
评分、协变量、Cox 设定与主模型（scripts/09 cohort_def + A6c）相同；评分在各自分析人群内标准化。
不并入任何合并分析，不替换任何已报告数值。"""
import os, sys, json, hashlib, importlib.util, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
warnings.filterwarnings("ignore")

sha = check_plan(); log("方案哈希校验通过", sha[:16])
DEV = f"{P}/docs/偏离与修订记录_外审回应A1-A6.md"; dev_txt = open(DEV, encoding="utf-8").read()
assert "## D5（2026-10-21" in dev_txt, "偏离记录中未找到 D5，拒绝运行"
dev_sha = hashlib.sha256(dev_txt.encode("utf-8")).hexdigest(); log("偏离记录 SHA-256", dev_sha[:16])
OUT = outdir("review_A6o")
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("l4", f"{P}/scripts/09_layer4_clinical.py"); l4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
MYE = "CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()
top = pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv", sep="\t").gene.tolist(); assert len(top) == 30

# 主模型（RFS）人群：用于界定 (b)/(c)，并复现表 S6a 的 GSE37745 主估计 0.78 (0.25–2.48)
d_rfs, g_rfs, *_ = l4.cohort_def("GSE37745")
assert len(d_rfs) == 34, len(d_rfs)

clin, expr = l4.read_series(f"{P}/raw/stageI_cohorts/GSE37745_series_matrix.txt.gz")
c = clin[(clin["histology"] == "adeno") & clin["tumor stage"].isin(["1a", "1b"])]
num = lambda x: pd.to_numeric(x, errors="coerce")
d = pd.DataFrame({"time": num(c["days to determined death status"]), "event": (c["dead"] == "yes").astype(int),
                  "age": num(c["age"]), "male": (c["gender"] == "male").astype(float),
                  "stage": (c["tumor stage"] == "1b").astype(float)}, index=c.index)
n_elig = len(d)
d = d.dropna(subset=["time"]); d = d[d.time > 0]
g = l4.to_gene(expr, "GPL570")
if g.stack().max() > 100: g = np.log2(g.clip(lower=1))
d = d.loc[[x for x in d.index if x in g.columns]]
log("分期合格", n_elig, "；OS 可用且有表达", len(d))
assert n_elig == 70, n_elig

def zr(m): return m.sub(m.mean(1), axis=0).div(m.std(1) + 1e-9, axis=0)
def score(gg, genes):
    pres = [x for x in genes if x in gg.index]; s = zr(gg.loc[pres].astype(float)).mean(0); return (s - s.mean()) / (s.std() + 1e-9)
def fit(dd, label):
    dd = dd.copy(); dd["age"] = dd.age.fillna(dd.age.median())
    gg = g[dd.index]; s = score(gg, top); m = score(gg, MYE)
    X0 = dd[["age", "male", "stage"]].copy(); X0["age"] = (X0.age - X0.age.mean()) / (X0.age.std() + 1e-9); X0 = X0.loc[:, X0.std() > 0]
    X = pd.concat([s.rename("score"), X0, m.rename("mye")], axis=1)
    r = PHReg(dd.time.values, X.values, status=dd.event.values, ties="breslow").fit(disp=False)
    b, se = float(r.params[0]), float(r.bse[0])
    out = dict(analysis=label, n=len(dd), events=int(dd.event.sum()), n_params=X.shape[1], covariates=",".join(X.columns),
               HR=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)), hi=float(np.exp(b + 1.96 * se)), beta=b, se=se,
               p_two_sided=float(2 * norm.sf(abs(b / se))), firth_needed=bool(dd.event.sum() < 20))
    log(label, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out.items()})
    return out

# 复现核对：RFS 主模型
dr = d_rfs.copy(); gg = g_rfs
s = score(gg, top); m = score(gg, MYE)
X0 = dr[["age", "male", "stage"]].copy(); X0["age"] = (X0.age - X0.age.mean()) / (X0.age.std() + 1e-9); X0 = X0.loc[:, X0.std() > 0]
rr = PHReg(dr.time.values, pd.concat([s, X0, m], axis=1).values, status=dr.event.values, ties="breslow").fit(disp=False)
rep = (round(float(np.exp(rr.params[0])), 2), round(float(np.exp(rr.params[0] - 1.96 * rr.bse[0])), 2), round(float(np.exp(rr.params[0] + 1.96 * rr.bse[0])), 2))
log("复现 RFS 主模型", rep, "（表 S6a：0.78, 0.25, 2.48）")

kept = set(d_rfs.index)
res = [fit(d, "(a) OS, all stage-eligible with OS data"),
       fit(d.loc[[x for x in d.index if x in kept]], "(b) OS, the 34 patients in the RFS analysis"),
       fit(d.loc[[x for x in d.index if x not in kept]], "(c) OS, the 36 patients excluded from the RFS analysis")]
pd.DataFrame(res).to_csv(OUT + "A6o_GSE37745_OS.tsv", sep="\t", index=False)
dump({"plan_sha256": sha, "deviation_record_sha256_at_run": dev_sha, "deviation": "D5 (post hoc)",
      "rfs_main_model_reproduced": rep, "rfs_main_model_reported": [0.78, 0.25, 2.48], "reproduces": rep == (0.78, 0.25, 2.48),
      "n_stage_eligible": n_elig, "results": res,
      "note": "Post hoc. Not pooled; does not replace any reported value. OS in early-stage LUAD includes non-cancer deaths."}, OUT + "A6o_results.json")
log("A6o 完成")
