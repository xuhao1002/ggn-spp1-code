#!/usr/bin/env python
"""第 5 层：可成药与长期安全（标准见 docs/项目进展与交接.md，2026-10-08 用户确认）。仅对上调基因判定；下调基因只报告 DGIdb 激动剂类药物（探索）。
可成药：DGIdb 有抑制类药物 — inhibitor/antagonist/blocker/antibody/negative modulator/inverse agonist/suppressor/antisense oligonucleotide（2026-10-08 修正：去掉 cleavage）；
        A = 已上市药物且来自 ChEMBL mechanism 或 TTD、作用类型为抑制；B = 其余（来自 ChEMBL/TTD/GuideToPharmacology/DTC）。
        类型缺失（NaN）或不明确（modulator/binder/other/unknown）且来自 ChEMBL/TTD/GuideToPharmacology/DTC 的相互作用计入，标「类型待核」。
        agonist/activator/positive modulator/potentiator/vaccine/immunotherapy 不计入（方向相反）。
安全：GTEx v8 心（左心室、心房）、肝、肾（皮质、髓质）、全血任一中位 TPM>=50 → 高表达否决；10–50 → 中等（标注，不否决）。
      HPA 25.1：心肌（cardiomyocytes）、肝细胞（hepatocytes）、肾（cells in glomeruli、cells in tubules 及各段肾小管）、骨髓（hematopoietic cells），
      可靠性 Approved/Supported/Enhanced 的记录中出现「High」→ 否决（含骨髓）；「Medium」→ 标注；无记录 → 标「蛋白证据缺失」，只看 GTEx。
      第 3 层的常见必需基因标「长期毒性风险」。
通过 = 可成药（A 或 B）且 GTEx 与 HPA 无否决。输出 results/10_layer5_safety_druggable/layer5.tsv、down_agonists.tsv"""
import os, zipfile, numpy as np, pandas as pd
P = os.path.expanduser("~/全部数据库/GGN_target"); S = f"{P}/raw/safety_druggability"; OUT = f"{P}/results/10_layer5_safety_druggable"
cand = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
l3 = pd.read_csv(f"{P}/results/08_layer3_depmap/layer3_depmap.tsv", sep="\t", index_col=0)

# ---- DGIdb ----
d = pd.read_csv(f"{S}/dgidb_interactions.tsv", sep="\t", low_memory=False)
d["itype"] = d["interaction_type"].fillna("").str.lower()
SRC = ["ChEMBL", "TTD", "GuideToPharmacology", "DTC"]
INHIB = {"inhibitor", "antagonist", "blocker", "antibody", "negative modulator", "inverse agonist", "suppressor", "antisense oligonucleotide"}   # 2026-10-08 修正：去掉 cleavage（降解底物，不是抑制该靶点）
UNSPEC = {"", "modulator", "binder", "other/unknown"}
AGON = {"agonist", "activator", "positive modulator", "potentiator"}
d = d[d["drug_name"].notna() | d["drug_claim_name"].notna()].copy()
d["drug"] = d["drug_name"].fillna(d["drug_claim_name"])
d["approved_b"] = d["approved"].astype(str).str.lower().eq("true")
inh = d[d["itype"].isin(INHIB)]
uns = d[d["itype"].isin(UNSPEC) & d["interaction_source_db_name"].isin(SRC)]
ag = d[d["itype"].isin(AGON)]

def summarize(g):
    """修正版（2026-10-08，用户批准）：不计 cleavage；A 级须有已上市药物且来自 ChEMBL mechanism 或 TTD 记录表明对该靶点的作用为抑制
    （即来源为 ChEMBL 或 TTD，且作用类型为抑制类）；其余有抑制类/类型待核记录（来自 ChEMBL/TTD/GuideToPharmacology/DTC）的为 B。"""
    a = inh[inh["gene_name"] == g]; u = uns[uns["gene_name"] == g]
    allr = pd.concat([a, u])
    if allr.empty: return dict(drug_grade="无", n_drugs=0, approved_drugs="", other_drugs="", type_to_verify=False)
    a_ok = a[a["approved_b"] & a["interaction_source_db_name"].isin(["ChEMBL", "TTD"])]
    appr = a_ok["drug"].str.upper().unique().tolist()
    other = allr[~allr["drug"].str.upper().isin(appr)]["drug"].str.upper().unique().tolist()
    grade = "A" if appr else "B"
    verify = a.empty                                   # 只有类型缺失/不明确的记录
    return dict(drug_grade=grade, n_drugs=len(set(appr) | set(other)), approved_drugs=";".join(appr[:6]), other_drugs=";".join(other[:6]),
                type_to_verify=bool(verify), explicit_inhibitory=not a.empty)

# ---- GTEx ----
g = pd.read_csv(f"{S}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
tis = ["Heart - Left Ventricle", "Heart - Atrial Appendage", "Liver", "Kidney - Cortex", "Kidney - Medulla", "Whole Blood"]
g = g.groupby("Description")[tis].max()
# ---- HPA ----
z = zipfile.ZipFile(f"{S}/normal_ihc_data.tsv.zip"); ih = pd.read_csv(z.open(z.namelist()[0]), sep="\t")
ih = ih[ih["Reliability"].isin(["Approved", "Supported", "Enhanced"])]
ct = {"Heart muscle": ["cardiomyocytes"], "Liver": ["hepatocytes"], "Bone marrow": ["hematopoietic cells"],
      "Kidney": ["cells in glomeruli", "cells in tubules", "proximal tubules (cell body)", "proximal tubules (microvilli)", "distal tubules", "collecting ducts"]}
ih = ih[ih.apply(lambda r: r["Cell type"] in ct.get(r["Tissue"], []), axis=1)]
hp = ih.groupby("Gene name")["Level"].apply(lambda x: "High" if (x == "High").any() else ("Medium" if (x == "Medium").any() else "Low/ND")).to_dict()
hp_cells = ih[ih["Level"] == "High"].groupby("Gene name").apply(lambda x: ";".join(sorted(set(x["Tissue"])))).to_dict()

rows = []
for gene, r in cand.iterrows():
    row = dict(gene=gene, direction=r["direction"], set=r["set"], in_main352=r["in_main352"], tier=r["tier"])
    if r["direction"] == "down":
        a = ag[ag["gene_name"] == gene]
        row.update(status="不判定（下调）", agonist_drugs=";".join(a["drug"].str.upper().unique()[:8]),
                   n_agonist_drugs=a["drug"].nunique(), agonist_approved=bool(a["approved_b"].any()))
        rows.append(row); continue
    row.update(summarize(gene))
    if gene in g.index:
        v = g.loc[gene]; mx = float(v.max()); row.update(gtex_max_tpm=round(mx, 1), gtex_max_tissue=v.idxmax(),
                      gtex_level="高表达（否决）" if mx >= 50 else ("中等" if mx >= 10 else "低"))
    else:
        row.update(gtex_max_tpm=np.nan, gtex_max_tissue="", gtex_level="无数据")
    h = hp.get(gene)
    row["hpa_level"] = {"High": "High（否决）", "Medium": "Medium（标注）", "Low/ND": "低/未检出"}.get(h, "蛋白证据缺失")
    row["hpa_high_tissues"] = hp_cells.get(gene, "")
    row["long_term_toxicity_risk"] = bool(l3["common_essential"].get(gene, False)) if "common_essential" in l3 else False
    druggable = row["drug_grade"] in ("A", "B")
    veto = row["gtex_level"].startswith("高") or row["hpa_level"].startswith("High")
    row["druggable"] = druggable; row["safety_veto"] = veto
    row["status"] = "通过" if (druggable and not veto) else "未通过"
    row["reason"] = "；".join([x for x, c in [("无可用抑制类药物", not druggable), ("GTEx 高表达", row["gtex_level"].startswith("高")), ("HPA High", row["hpa_level"].startswith("High"))] if c])
    rows.append(row)
out = pd.DataFrame(rows).set_index("gene")
out.to_csv(f"{OUT}/layer5.tsv", sep="\t")
out[out.direction == "down"].dropna(subset=["n_agonist_drugs"]).query("n_agonist_drugs>0").to_csv(f"{OUT}/down_agonists.tsv", sep="\t")
up = out[out.direction == "up"]
print("上调", len(up), up["status"].value_counts().to_dict())
for lab, sub in [("预设 352", up[up.in_main352]), ("683", up)]:
    print(lab, len(sub), sub["status"].value_counts().to_dict(), "；可成药 A", int((sub.drug_grade == "A").sum()), "B", int((sub.drug_grade == "B").sum()),
          "；类型待核", int(sub.type_to_verify.fillna(False).sum()), "；GTEx 否决", int(sub.gtex_level.str.startswith("高").sum()), "；HPA High", int(sub.hpa_level.str.startswith("High").sum()))
dn = out[out.direction == "down"]
print("下调", len(dn), "；有激动剂类药物的", int((dn["n_agonist_drugs"] > 0).sum()))
