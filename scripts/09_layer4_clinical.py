#!/usr/bin/env python
"""
第 4 层：I 期队列的临床相关（标准见 docs/项目进展与交接.md「第 3–5 层判定标准」，2026-10-08 用户确认）。
- 队列：GSE31210、GSE50081、GSE30219、GSE37745、GSE68465、GSE72094（探针按平台注释转基因名，多探针取平均表达最高者）、TCGA-LUAD（Ensembl→基因名用 GENCODE v44）。
- 入组：肺腺癌、病理 I 期；GSE31210/GSE50081 先筛 I 期（IA、IB）；GSE30219、GSE68465 无分期，用 pT1–2 且 pN0（GSE30219 另 M0）近似。
- 终点：RFS/DFS 优先（GSE31210、GSE50081、GSE30219、GSE37745、GSE68465），否则 OS（GSE72094、TCGA）。
- 模型：Cox（statsmodels PHReg，Breslow），基因表达在队列内 z 分数，协变量=年龄、性别、分期（IB/T2 对 IA/T1）；每队列给 log HR/se；单侧 P 按预期方向。
- 通过（上调基因）：>=2 个队列 HR>1 且单侧 P<0.05；随机效应（DerSimonian–Laird）meta 单侧 P 的 BH q<0.05（在上调基因中校正）；无队列显著反向（HR<1 且单侧 P<0.05）。
- 下调基因只报告（预期 HR<1），不判定。
输出 results/09_layer4_clinical/：layer4_clinical.tsv、cohort_summary.tsv
"""
import os, re, gzip, io, warnings
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/09_layer4_clinical"; os.makedirs(OUT, exist_ok=True)
cand = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
GENES = list(cand.index)

def read_series(f):
    """按样本逐个解析「键: 值」（同一行里不同样本的键可能不同）；同一样本内重复的键加 #1、#2 后缀。"""
    ids, rows = None, []
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"): ids = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_characteristics_ch1"):
                rows.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    recs = []
    for j in range(len(ids)):
        d = {}
        for r in rows:
            if j >= len(r) or ":" not in r[j]: continue
            k, v = r[j].split(":", 1); k, v = k.strip(), v.strip()
            kk, i = k, 1
            while kk in d: kk = f"{k}#{i}"; i += 1
            d[kk] = v
        recs.append(d)
    clin = pd.DataFrame(recs, index=ids)
    expr = pd.read_csv(f, sep="\t", comment="!", index_col=0, compression="gzip")
    expr.columns = [c.strip('"') for c in expr.columns]
    return clin, expr

def probe_map(plat):
    if plat == "GPL15048":
        with gzip.open(f"{P}/raw/annotation/platforms/GPL15048.txt.gz", "rt") as fh:
            txt = fh.read()
        t = txt.split("!platform_table_begin\n")[1].split("\n!platform_table_end")[0]
        d = pd.read_csv(io.StringIO(t), sep="\t", dtype=str); return d.set_index("ID")["GeneSymbol"]
    with gzip.open(f"{P}/raw/annotation/platforms/{plat}.annot.gz", "rt", errors="ignore") as fh:
        lines = fh.read().split("\n")
    i = [k for k, l in enumerate(lines) if l.startswith("!platform_table_begin")][0]
    j = [k for k, l in enumerate(lines) if l.startswith("!platform_table_end")]
    t = "\n".join(lines[i + 1:(j[0] if j else len(lines))])
    d = pd.read_csv(io.StringIO(t), sep="\t", dtype=str)
    return d.set_index("ID")["Gene symbol"]

def to_gene(expr, plat):
    pm = probe_map(plat)
    s = pm.reindex(expr.index)
    ok = s.notna() & (s != "") & ~s.astype(str).str.contains("///")
    e = expr[ok.values].copy(); e["gene"] = s[ok].values
    e["mean"] = e.drop(columns="gene").mean(axis=1)
    e = e.sort_values("mean", ascending=False).drop_duplicates("gene").set_index("gene").drop(columns="mean")
    return e

def num(x): return pd.to_numeric(x, errors="coerce")

def cohort_def(name):
    f = f"{P}/raw/stageI_cohorts/{name}_series_matrix.txt.gz"
    clin, expr = read_series(f)
    if name == "GSE31210":
        c = clin[clin["pathological stage"].isin(["IA", "IB"])]
        d = pd.DataFrame({"time": num(c["months before relapse/censor"]), "event": (c["relapse"] == "relapsed").astype(int),
                          "age": num(c["age (years)"]), "male": (c["gender"] == "male").astype(float),
                          "stage": (c["pathological stage"] == "IB").astype(float)}, index=c.index)
        plat, endpoint = "GPL570", "RFS"
    elif name == "GSE50081":
        k = [x for x in clin.columns if x.startswith("recurrence")][0]
        c = clin[(clin["histology"] == "adenocarcinoma") & clin["Stage"].isin(["1A", "1B"])]
        d = pd.DataFrame({"time": num(c["disease-free survival time"]), "event": (c[k] == "Y").astype(int),
                          "age": num(c["age"]), "male": (c["Sex"] == "M").astype(float),
                          "stage": (c["Stage"] == "1B").astype(float)}, index=c.index)
        plat, endpoint = "GPL570", "DFS"
    elif name == "GSE72094":
        c = clin[clin["Stage"].astype(str).isin(["1A", "1B"])]
        d = pd.DataFrame({"time": num(c["survival_time_in_days"]), "event": (c["vital_status"] == "Dead").astype(int),
                          "age": num(c["age_at_diagnosis"]), "male": (c["gender"] == "M").astype(float),
                          "stage": (c["Stage"] == "1B").astype(float)}, index=c.index)
        plat, endpoint = "GPL15048", "OS"
    elif name == "GSE37745":
        c = clin[(clin["histology"] == "adeno") & clin["tumor stage"].isin(["1a", "1b"])]
        d = pd.DataFrame({"time": num(c["days to recurrence / to last visit"]), "event": (c["recurrence"] == "yes").astype(int),
                          "age": num(c["age"]), "male": (c["gender"] == "male").astype(float),
                          "stage": (c["tumor stage"] == "1b").astype(float)}, index=c.index)
        plat, endpoint = "GPL570", "RFS"
    elif name == "GSE30219":
        k = [x for x in clin.columns if x.startswith("relapse")][0]
        c = clin[(clin["histology"] == "ADC") & clin["pt stage"].isin(["T1", "T2"]) & (clin["pn stage"] == "N0") & (clin["pm stage"] == "M0")]
        d = pd.DataFrame({"time": num(c["disease free survival in months"]), "event": num(c[k]),
                          "age": num(c["age at surgery"]), "male": (c["gender"] == "M").astype(float),
                          "stage": (c["pt stage"] == "T2").astype(float)}, index=c.index)
        plat, endpoint = "GPL570", "DFS"
    elif name == "GSE68465":
        c = clin[clin["disease_stage"].isin(["pN0pT1", "pN0pT2"])]
        prog = c["first_progression_or_relapse"] == "Yes"
        t = np.where(prog, num(c["months_to_first_progression"]), num(c["months_to_last_contact_or_death"]))
        d = pd.DataFrame({"time": t, "event": prog.astype(int), "age": num(c["age"]), "male": (c["Sex"] == "Male").astype(float),
                          "stage": (c["disease_stage"] == "pN0pT2").astype(float)}, index=c.index)
        plat, endpoint = "GPL96", "RFS"
    g = to_gene(expr, plat)
    d = d.dropna(subset=["time", "event"]); d = d[d["time"] > 0]
    d["age"] = d["age"].fillna(d["age"].median())
    raw_scale = g.stack().max() > 100          # MAS5 原始强度 → log2
    if raw_scale: g = np.log2(g.clip(lower=1))
    g = g[[x for x in d.index if x in g.columns]]
    return d.loc[g.columns], g, plat, endpoint, raw_scale

def tcga():
    clin = pd.read_csv(f"{P}/raw/TCGA_LUAD/TCGA-LUAD.clinical.tsv.gz", sep="\t", index_col=0)
    sv = pd.read_csv(f"{P}/raw/TCGA_LUAD/TCGA-LUAD.survival.tsv.gz", sep="\t", index_col=0)
    ed = [c for c in clin.columns if "ajcc_staging_system_edition" in c]
    c = clin[clin["ajcc_pathologic_stage.diagnoses"].isin(["Stage I", "Stage IA", "Stage IB"]) & (clin["sample_type.samples"] == "Primary Tumor")]
    c = c[c.index.isin(sv.index)]
    d = pd.DataFrame({"time": sv.loc[c.index, "OS.time"].astype(float), "event": sv.loc[c.index, "OS"].astype(int),
                      "age": num(c["age_at_index.demographic"]), "male": (c["gender.demographic"] == "male").astype(float),
                      "stage": (c["ajcc_pathologic_stage.diagnoses"] == "Stage IB").astype(float)}, index=c.index)
    d = d.dropna(subset=["time"]); d = d[d["time"] > 0]; d["age"] = d["age"].fillna(d["age"].median())
    if ed: print("TCGA AJCC 版本分布：", c[ed[0]].value_counts().to_dict())
    # 基因表达
    gtf = pd.read_csv(f"{P}/raw/annotation/gencode.v44.basic.annotation.gtf.gz", sep="\t", comment="#", header=None, usecols=[2, 8], names=["f", "a"], low_memory=False)
    gtf = gtf[gtf["f"] == "gene"]
    mp = pd.Series(gtf["a"].str.extract(r'gene_name "([^"]+)"')[0].values, index=gtf["a"].str.extract(r'gene_id "([^".]+)')[0].values)
    mp = mp[~mp.index.duplicated()]
    hdr = pd.read_csv(f"{P}/raw/TCGA_LUAD/TCGA-LUAD.star_tpm.tsv.gz", sep="\t", nrows=0, compression="gzip").columns
    cols = [h for h in hdr[1:] if h in d.index]
    e = pd.read_csv(f"{P}/raw/TCGA_LUAD/TCGA-LUAD.star_tpm.tsv.gz", sep="\t", usecols=["Ensembl_ID"] + cols, index_col=0, compression="gzip")
    e.index = e.index.str.replace(r"\.\d+$", "", regex=True)
    e["gene"] = mp.reindex(e.index).values; e = e[e["gene"].notna()]
    e["mean"] = e.drop(columns="gene").mean(axis=1)
    e = e.sort_values("mean", ascending=False).drop_duplicates("gene").set_index("gene").drop(columns="mean")
    return d.loc[cols], e, "RNA-seq", "OS", False

def cox_all(d, g, genes):
    X0 = d[["age", "male", "stage"]].copy()
    X0["age"] = (X0["age"] - X0["age"].mean()) / (X0["age"].std() + 1e-9)
    X0 = X0.loc[:, X0.std() > 0]
    out = {}
    for gene in genes:
        if gene not in g.index: continue
        x = g.loc[gene].astype(float)
        if x.std() == 0 or x.isna().any(): continue
        z = (x - x.mean()) / x.std()
        X = pd.concat([z.rename("g"), X0], axis=1).values
        try:
            r = PHReg(d["time"].values, X, status=d["event"].values, ties="breslow").fit(disp=False)
            out[gene] = (float(r.params[0]), float(r.bse[0]))
        except Exception:
            continue
    return out

def main():
    res, summ = {}, []
    for name in ["GSE31210", "GSE50081", "GSE30219", "GSE37745", "GSE68465", "GSE72094", "TCGA"]:
        d, g, plat, endpoint, rawscale = tcga() if name == "TCGA" else cohort_def(name)
        r = cox_all(d, g, GENES)
        res[name] = r
        summ.append(dict(cohort=name, platform=plat, endpoint=endpoint, n=len(d), events=int(d["event"].sum()),
                         n_IB_or_T2=int(d["stage"].sum()), log2_transformed=rawscale, genes_tested=len(r)))
        print(summ[-1], flush=True)
    pd.DataFrame(summ).to_csv(f"{OUT}/cohort_summary.tsv", sep="\t", index=False)
    rows = []
    for gene, row in cand.iterrows():
        up = row["direction"] == "up"
        d = dict(gene=gene, direction=row["direction"], set=row["set"], in_main352=row["in_main352"], tier=row["tier"])
        lh, se = [], []
        n_sup = n_opp = 0
        for name, r in res.items():
            if gene not in r: d[f"{name}_logHR"] = np.nan; continue
            b, s = r[gene]; z = b / s
            p_exp = (1 - norm.cdf(z)) if up else norm.cdf(z)          # 预期方向的单侧 P
            p_opp = norm.cdf(z) if up else (1 - norm.cdf(z))
            d[f"{name}_logHR"] = round(b, 4); d[f"{name}_p_exp"] = round(p_exp, 4)
            lh.append(b); se.append(s)
            n_sup += p_exp < 0.05; n_opp += p_opp < 0.05
        d["n_cohorts"] = len(lh); d["n_support"] = n_sup; d["n_opposite_sig"] = n_opp
        if len(lh) >= 2:
            lh, se = np.array(lh), np.array(se); w = 1 / se**2
            fe = (w * lh).sum() / w.sum(); Q = (w * (lh - fe) ** 2).sum()
            tau2 = max(0, (Q - (len(lh) - 1)) / (w.sum() - (w**2).sum() / w.sum()))
            w2 = 1 / (se**2 + tau2); mu = (w2 * lh).sum() / w2.sum(); se_mu = np.sqrt(1 / w2.sum())
            z = mu / se_mu
            d["meta_logHR"] = round(mu, 4); d["meta_HR"] = round(float(np.exp(mu)), 3)
            d["meta_p_exp"] = (1 - norm.cdf(z)) if up else norm.cdf(z); d["I2"] = round(max(0, (Q - (len(lh) - 1)) / Q) if Q > 0 else 0, 3)
        rows.append(d)
    out = pd.DataFrame(rows).set_index("gene")
    mask = (out["direction"] == "up") & out["meta_p_exp"].notna()
    out["meta_q"] = np.nan
    out.loc[mask, "meta_q"] = multipletests(out.loc[mask, "meta_p_exp"], method="fdr_bh")[1]
    ok = (out["direction"] == "up") & (out["n_support"] >= 2) & (out["meta_q"] < 0.05) & (out["n_opposite_sig"] == 0)
    out["status"] = np.where(out["direction"] == "down", "不判定（下调，仅报告）",
                     np.where(out["n_cohorts"] < 2, "无法检验", np.where(ok, "通过", "未通过")))
    out.to_csv(f"{OUT}/layer4_clinical.tsv", sep="\t")
    up = out[out.direction == "up"]
    print("上调基因", len(up), up["status"].value_counts().to_dict())
    for lab, sub in [("预设 352", up[up.in_main352]), ("683", up)]:
        print(lab, len(sub), sub["status"].value_counts().to_dict())
    p = up[up.status == "通过"].sort_values("meta_q")
    print("通过者（按 meta q）：", ", ".join(f"{g}(HR{r.meta_HR},q{r.meta_q:.3f},{int(r.n_support)}/{int(r.n_cohorts)})" for g, r in p.head(40).iterrows()))
    dn = out[out.direction == "down"]
    print("下调基因：", len(dn), "；>=2 队列 HR<1 且单侧 P<0.05：", int((dn["n_support"] >= 2).sum()))

if __name__ == "__main__":
    main()
