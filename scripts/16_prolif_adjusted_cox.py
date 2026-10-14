#!/usr/bin/env python
"""探索性（不参与判定）：第 4 层 Cox 加增殖评分协变量（标准见 docs/项目进展与交接.md，2026-10-08 用户确认）。
增殖评分 = 队列内 Tirosh 96 个基因（S 43 + G2/M 53）z 分数的平均（队列中缺失的基因忽略）；被检验基因若属于该基因集，则不调整，标「属于增殖基因集」。
报告调整前后 meta HR 与 q；「独立于增殖」=调整后仍满足第 4 层通过标准（>=2 队列 HR>1 且单侧 P<0.05；meta q<0.05；无队列显著反向）。
用法：python scripts/16_prolif_adjusted_cox.py → results/16_prolif_adjusted/"""
import os, sys, importlib.util
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.expanduser("~/全部数据库/GGN_target/scripts"))
import prolif_genes
spec = importlib.util.spec_from_file_location("l4", os.path.expanduser("~/全部数据库/GGN_target/scripts/09_layer4_clinical.py"))
l4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
P = l4.P; OUT = f"{P}/results/16_prolif_adjusted"; os.makedirs(OUT, exist_ok=True)
GENES = l4.GENES; TIROSH = set(prolif_genes.TIROSH)

def cox_adj(d, g, genes):
    present = [x for x in prolif_genes.TIROSH if x in g.index]
    z = g.loc[present].astype(float); z = (z.sub(z.mean(axis=1), axis=0)).div(z.std(axis=1) + 1e-9, axis=0)
    score = z.mean(axis=0).reindex(d.index); score = (score - score.mean()) / (score.std() + 1e-9)
    X0 = d[["age", "male", "stage"]].copy(); X0["age"] = (X0["age"] - X0["age"].mean()) / (X0["age"].std() + 1e-9)
    X0 = X0.loc[:, X0.std() > 0]; X0["prolif"] = score
    out = {}
    for gene in genes:
        if gene not in g.index or gene in TIROSH: continue
        x = g.loc[gene].astype(float)
        if x.std() == 0 or x.isna().any(): continue
        zz = (x - x.mean()) / x.std()
        X = pd.concat([zz.rename("g"), X0], axis=1).values
        try:
            r = PHReg(d["time"].values, X, status=d["event"].values, ties="breslow").fit(disp=False)
            out[gene] = (float(r.params[0]), float(r.bse[0]))
        except Exception: continue
    return out, len(present)

def main():
    cand = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
    base = pd.read_csv(f"{P}/results/09_layer4_clinical/layer4_clinical.tsv", sep="\t", index_col=0)
    res = {}
    for name in ["GSE31210", "GSE50081", "GSE30219", "GSE37745", "GSE68465", "GSE72094", "TCGA"]:
        d, g, plat, endpoint, rs = l4.tcga() if name == "TCGA" else l4.cohort_def(name)
        r, npres = cox_adj(d, g, GENES); res[name] = r
        print(name, "增殖评分基因", npres, "；检验基因", len(r), flush=True)
    rows = []
    for gene, row in cand.iterrows():
        up = row["direction"] == "up"; d = dict(gene=gene, direction=row["direction"], set=row["set"], in_main352=row["in_main352"], in_prolif_set=gene in TIROSH)
        lh, se = [], []; n_sup = n_opp = 0
        for name, r in res.items():
            if gene not in r: continue
            b, s = r[gene]; z = b / s
            p_exp = (1 - norm.cdf(z)) if up else norm.cdf(z); p_opp = norm.cdf(z) if up else (1 - norm.cdf(z))
            lh.append(b); se.append(s); n_sup += p_exp < 0.05; n_opp += p_opp < 0.05
        d["n_cohorts_adj"] = len(lh); d["n_support_adj"] = n_sup; d["n_opposite_adj"] = n_opp
        if len(lh) >= 2:
            lh, se = np.array(lh), np.array(se); w = 1 / se**2; fe = (w * lh).sum() / w.sum(); Q = (w * (lh - fe) ** 2).sum()
            tau2 = max(0, (Q - (len(lh) - 1)) / (w.sum() - (w**2).sum() / w.sum())); w2 = 1 / (se**2 + tau2)
            mu = (w2 * lh).sum() / w2.sum(); z = mu / np.sqrt(1 / w2.sum())
            d["meta_HR_adj"] = round(float(np.exp(mu)), 3); d["meta_p_adj"] = (1 - norm.cdf(z)) if up else norm.cdf(z)
        rows.append(d)
    out = pd.DataFrame(rows).set_index("gene")
    m = (out["direction"] == "up") & out["meta_p_adj"].notna()
    out["meta_q_adj"] = np.nan; out.loc[m, "meta_q_adj"] = multipletests(out.loc[m, "meta_p_adj"], method="fdr_bh")[1]
    ok = (out["direction"] == "up") & (out["n_support_adj"] >= 2) & (out["meta_q_adj"] < 0.05) & (out["n_opposite_adj"] == 0)
    out["status_adj"] = np.where(out["direction"] == "down", "不判定", np.where(out["in_prolif_set"], "属于增殖基因集", np.where(out["n_cohorts_adj"] < 2, "无法检验", np.where(ok, "通过", "未通过"))))
    out["status_orig"] = base["status"].reindex(out.index); out["meta_HR_orig"] = base["meta_HR"]; out["meta_q_orig"] = base["meta_q"]
    out.to_csv(f"{OUT}/layer4_prolif_adjusted.tsv", sep="\t")
    up = out[out.direction == "up"]
    for lab, sub in [("预设 352", up[up.in_main352]), ("683", up)]:
        orig = sub[sub.status_orig == "通过"]
        print(lab, f"上调 {len(sub)}；调整前通过 {len(orig)}；调整后通过 {int((sub.status_adj=='通过').sum())}；原通过者中仍通过 {int((orig.status_adj=='通过').sum())}，属于增殖基因集 {int((orig.status_adj=='属于增殖基因集').sum())}，不再通过 {int((orig.status_adj=='未通过').sum())}")
    keep = up[(up.status_adj == "通过") & (up.status_orig == "通过")].sort_values("meta_q_adj")
    print("独立于增殖（调整前后都通过）：", ", ".join(f"{g}(HR{r.meta_HR_orig}→{r.meta_HR_adj},q{r.meta_q_adj:.3f})" for g, r in keep.head(40).iterrows()))
    newp = up[(up.status_adj == "通过") & (up.status_orig != "通过")]
    print("调整后新通过（调整前未通过）：", ", ".join(newp.index[:30]))

if __name__ == "__main__":
    main()
