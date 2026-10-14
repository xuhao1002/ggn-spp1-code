#!/usr/bin/env python
"""
微环境第 4、5 层（标准见 docs/项目进展与交接.md「补充分析判定标准」，2026-10-08 用户确认）。仅对通过第 1 层外部验证的上调基因（阻断靶点候选）判定；下调基因不判定。
第 3 层不适用，按设计无法满足「通过全部五层」：报告为「通过第 1、4、5 层 + 文献证据」，不计入成功标准；文献证据人工核查。
- 第 4 层：沿用 7 个 I 期队列与 Cox（年龄、性别、分期），额外加该细胞群的标志基因平均 z 分数为协变量（属于该群标志基因的基因不调整，标注）；
  通过=>=2 队列 HR>1 且单侧 P<0.05、meta q<0.05（BH，在该细胞群内检验的基因中校正）、无队列显著反向。
- 第 5 层：可成药（DGIdb 抑制类或类型待核；A/B 分级同上皮）且安全：心（左心室、心房）、肝、肾 GTEx>=50 TPM 或 HPA High（心肌、肝细胞、肾各段）否决；
  非免疫群（成纤维、内皮、平滑肌）另加全血>=50 TPM 与 HPA 骨髓造血细胞 High 否决；
  免疫群（T/NK、髓系、B/浆）的血液/骨髓否决改为：目标群表达 >= 其他免疫谱系最高者的 3 倍（Peng 样本级中位 log1p 表达，同来源过滤表）且 HPA 骨髓造血细胞无 High。
输出 results/15_microenv_l45/
"""
import os, sys, importlib.util, zipfile
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/15_microenv_l45"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("l4", f"{P}/scripts/09_layer4_clinical.py"); l4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
MARK = {"T/NK": "CD3E CD3D TRAC CD2 NKG7 GNLY KLRD1 KLRF1".split(),
        "Myeloid": "CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split(),
        "B/Plasma": "MS4A1 CD79A BANK1 JCHAIN MZB1 IGHG1".split(),
        "Fibroblast": "COL1A1 COL1A2 DCN LUM".split(), "Endothelial": "PECAM1 VWF CLDN5 CDH5".split(), "SMC/pericyte": "ACTA2 MYH11 RGS5 PDGFRB".split()}
IMMUNE = ["T/NK", "Myeloid", "B/Plasma"]
cand = pd.read_csv(f"{P}/results/14_microenv/microenv_layer1_candidates.tsv", sep="\t", index_col=0)
cand = cand[cand["L1_pass"] & (cand["log2FoldChange"] > 0)].copy()
cand["gene"] = cand.index
print("通过第 1 层的上调（基因×细胞群）", len(cand), cand["lineage"].value_counts().to_dict())
genes_all = sorted(set(cand.index))

# ---------- 第 4 层 ----------
cohorts = {}
for name in ["GSE31210", "GSE50081", "GSE30219", "GSE37745", "GSE68465", "GSE72094", "TCGA"]:
    cohorts[name] = l4.tcga() if name == "TCGA" else l4.cohort_def(name)
    print("队列", name, len(cohorts[name][0]), flush=True)

def cox_lineage(d, g, genes, markers):
    pres = [x for x in markers if x in g.index]
    z = g.loc[pres].astype(float); z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
    sc_ = z.mean(axis=0).reindex(d.index); sc_ = (sc_ - sc_.mean()) / (sc_.std() + 1e-9)
    X0 = d[["age", "male", "stage"]].copy(); X0["age"] = (X0["age"] - X0["age"].mean()) / (X0["age"].std() + 1e-9)
    X0 = X0.loc[:, X0.std() > 0]; X0["lin"] = sc_
    out = {}
    for gene in genes:
        if gene not in g.index or gene in markers: continue
        x = g.loc[gene].astype(float)
        if x.std() == 0 or x.isna().any(): continue
        X = pd.concat([((x - x.mean()) / x.std()).rename("g"), X0], axis=1).values
        try:
            r = PHReg(d["time"].values, X, status=d["event"].values, ties="breslow").fit(disp=False); out[gene] = (float(r.params[0]), float(r.bse[0]))
        except Exception: pass
    return out
rows = []
for L in MARK:
    gl = list(cand[cand.lineage == L].index)
    if not gl: continue
    res = {n: cox_lineage(c[0], c[1], gl, MARK[L]) for n, c in cohorts.items()}
    for gene in gl:
        d = dict(gene=gene, lineage=L, in_marker_set=gene in MARK[L]); lh, se = [], []; n_sup = n_opp = 0
        for n, r in res.items():
            if gene not in r: continue
            b, s = r[gene]; z = b / s; lh.append(b); se.append(s); n_sup += (1 - norm.cdf(z)) < 0.05; n_opp += norm.cdf(z) < 0.05
        d.update(n_cohorts=len(lh), n_support=n_sup, n_opposite_sig=n_opp)
        if len(lh) >= 2:
            lh, se = np.array(lh), np.array(se); w = 1 / se**2; fe = (w * lh).sum() / w.sum(); Q = (w * (lh - fe) ** 2).sum()
            tau2 = max(0, (Q - (len(lh) - 1)) / (w.sum() - (w**2).sum() / w.sum())); w2 = 1 / (se**2 + tau2); mu = (w2 * lh).sum() / w2.sum()
            d["meta_HR"] = round(float(np.exp(mu)), 3); d["meta_p"] = 1 - norm.cdf(mu / np.sqrt(1 / w2.sum()))
        rows.append(d)
l4t = pd.DataFrame(rows)
l4t["meta_q"] = np.nan
for L in MARK:
    m = (l4t.lineage == L) & l4t.meta_p.notna()
    if m.sum(): l4t.loc[m, "meta_q"] = multipletests(l4t.loc[m, "meta_p"], method="fdr_bh")[1]
l4t["L4_pass"] = (l4t["n_support"] >= 2) & (l4t["meta_q"] < 0.05) & (l4t["n_opposite_sig"] == 0) & ~l4t["in_marker_set"]

# ---------- 第 5 层 ----------
S = f"{P}/raw/safety_druggability"
d = pd.read_csv(f"{S}/dgidb_interactions.tsv", sep="\t", low_memory=False); d["itype"] = d["interaction_type"].fillna("").str.lower()
SRC = ["ChEMBL", "TTD", "GuideToPharmacology", "DTC"]
INHIB = {"inhibitor", "antagonist", "blocker", "antibody", "negative modulator", "inverse agonist", "suppressor", "antisense oligonucleotide"}   # 2026-10-08 修正：去掉 cleavage
UNSPEC = {"", "modulator", "binder", "other/unknown"}
d = d[d["drug_name"].notna() | d["drug_claim_name"].notna()].copy(); d["drug"] = d["drug_name"].fillna(d["drug_claim_name"]); d["appr"] = d["approved"].astype(str).str.lower().eq("true")
inh = d[d["itype"].isin(INHIB)]; uns = d[d["itype"].isin(UNSPEC) & d["interaction_source_db_name"].isin(SRC)]
g = pd.read_csv(f"{S}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
core = ["Heart - Left Ventricle", "Heart - Atrial Appendage", "Liver", "Kidney - Cortex", "Kidney - Medulla"]
g = g.groupby("Description")[core + ["Whole Blood"]].max()
z = zipfile.ZipFile(f"{S}/normal_ihc_data.tsv.zip"); ih = pd.read_csv(z.open(z.namelist()[0]), sep="\t"); ih = ih[ih["Reliability"].isin(["Approved", "Supported", "Enhanced"])]
ctc = {"Heart muscle": ["cardiomyocytes"], "Liver": ["hepatocytes"], "Kidney": ["cells in glomeruli", "cells in tubules", "proximal tubules (cell body)", "proximal tubules (microvilli)", "distal tubules", "collecting ducts"]}
def hpa_high(tissues_ct):
    x = ih[ih.apply(lambda r: r["Cell type"] in tissues_ct.get(r["Tissue"], []), axis=1)]
    return set(x.loc[x["Level"] == "High", "Gene name"])
hp_core = hpa_high(ctc); hp_bm = hpa_high({"Bone marrow": ["hematopoietic cells"]})
ratio_p = pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv", sep="\t", index_col=0)
imm_cols = {"T/NK": "med_T/NK", "Myeloid": "med_Myeloid", "B/Plasma": "med_B/Plasma"}; other_imm = ["med_T/NK", "med_Myeloid", "med_B/Plasma", "med_Mast"]
rows5 = []
for _, r in l4t.iterrows():
    gene, L = r["gene"], r["lineage"]; row = dict(gene=gene, lineage=L)
    a = inh[inh["gene_name"] == gene]; u = uns[uns["gene_name"] == gene]; allr = pd.concat([a, u])
    a_ok = a[a["appr"] & a["interaction_source_db_name"].isin(["ChEMBL", "TTD"])]          # 2026-10-08 修正：A 级须来自 ChEMBL/TTD 的抑制类记录
    appr = a_ok["drug"].str.upper().unique().tolist(); oth = allr[~allr["drug"].str.upper().isin(appr)]["drug"].str.upper().unique().tolist()
    row["drug_grade"] = "无" if allr.empty else ("A" if appr else "B"); row["drugs"] = ";".join((appr + oth)[:6])
    row["type_to_verify"] = (not allr.empty) and a.empty
    tpm = g.loc[gene] if gene in g.index else None
    row["gtex_core_max"] = float(tpm[core].max()) if tpm is not None else np.nan; row["gtex_blood"] = float(tpm["Whole Blood"]) if tpm is not None else np.nan
    veto_core = (row["gtex_core_max"] >= 50) or (gene in hp_core)
    if L in IMMUNE:
        others = [c for c in other_imm if c != imm_cols[L]]
        mL = ratio_p.loc[gene, imm_cols[L]] if gene in ratio_p.index else np.nan; mo = ratio_p.loc[gene, others].max() if gene in ratio_p.index else np.nan
        row["lineage_fold_vs_other_immune"] = round(mL / mo, 2) if (mo and mo > 0) else (np.inf if mL and mL > 0 else np.nan)
        bm_rule = (row["lineage_fold_vs_other_immune"] >= 3) and (gene not in hp_bm)
        row["blood_marrow_rule"] = bool(bm_rule)
        veto = veto_core or not bm_rule
    else:
        bm_veto = (row["gtex_blood"] >= 50) or (gene in hp_bm); row["blood_marrow_rule"] = not bm_veto
        veto = veto_core or bm_veto
    row["safety_veto"] = bool(veto); row["L5_pass"] = (row["drug_grade"] in ("A", "B")) and not veto
    rows5.append(row)
l5t = pd.DataFrame(rows5)
out = l4t.merge(l5t, on=["gene", "lineage"]).merge(cand[["lineage", "log2FoldChange", "padj", "GSE164789_pass", "GSE189357_pass"]].reset_index(drop=True).assign(gene=cand.index.values), on=["gene", "lineage"], how="left")
out["pass_1_4_5"] = out["L4_pass"] & out["L5_pass"]
out.to_csv(f"{OUT}/microenv_layers45.tsv", sep="\t", index=False)
print(out.groupby("lineage").agg(L1=("gene", "size"), L4=("L4_pass", "sum"), L5=("L5_pass", "sum"), both=("pass_1_4_5", "sum")))
print("通过第 1、4、5 层（基因×群）：", int(out.pass_1_4_5.sum()), "；基因数", out[out.pass_1_4_5].gene.nunique())
for L, sub in out[out.pass_1_4_5].groupby("lineage"):
    print(L, ", ".join(f"{r.gene}(LFC{r.log2FoldChange:+.1f},HR{r.meta_HR},{r.drug_grade})" for _, r in sub.sort_values("meta_q").iterrows()))
