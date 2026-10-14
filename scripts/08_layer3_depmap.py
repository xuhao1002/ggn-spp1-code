#!/usr/bin/env python
"""第 3 层：DepMap（方案 A，标准见 docs/项目进展与交接.md「第 3–5 层判定标准」）。仅对上调基因判定；下调基因不判定。
- 依赖细胞系：Chronos<=-0.5；常见必需：>=90% 全部 CRISPR 模型依赖；肺腺癌 = Model.csv 中 OncotreeSubtype=='Lung Adenocarcinoma'
- 通过：肺腺癌系>=20% 依赖、非常见必需、肺腺癌系表达（log2(TPM+1)>=1 的系>=50%）
- 否决证据：常见必需；敲除促进生长（Chronos>=+0.3 的肺腺癌系>=20%）
- 无存活依赖（表达但依赖<20%）= 未通过，入「待补充池」（后续补 LINCS L1000 敲低特征）；基因不在 CRISPR 矩阵 = 无法检验
输出 results/08_layer3_depmap/layer3_depmap.tsv"""
import os, re, numpy as np, pandas as pd
B = os.path.expanduser("~/全部数据库"); P = f"{B}/GGN_target"
cand = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
m = pd.read_csv(f"{B}/GGN_target/raw/DepMap/Model.csv", low_memory=False)
luad_ids = set(m.loc[m["OncotreeSubtype"] == "Lung Adenocarcinoma", "ModelID"])
cr = pd.read_csv(f"{P}/raw/DepMap/CRISPRGeneEffect.csv", index_col=0)
cr.columns = [re.sub(r" \(\d+\)$", "", c) for c in cr.columns]
ex = pd.read_csv(f"{P}/raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", index_col=0, low_memory=False)
ex.columns = [re.sub(r" \(\d+\)$", "", c) for c in ex.columns]
if "ModelID" in ex.columns: ex = ex.set_index("ModelID")
l_cr = cr.index.intersection(luad_ids); l_ex = ex.index.intersection(luad_ids)
print("CRISPR 模型", len(cr), "；其中肺腺癌", len(l_cr), "；表达数据肺腺癌", len(l_ex))
rows = []
for g, r in cand.iterrows():
    d = dict(gene=g, direction=r["direction"], set=r["set"], in_main352=r["in_main352"], tier=r["tier"])
    if r["direction"] == "down":
        d.update(status="不判定（下调，入缺失/恢复清单）", reason=""); rows.append(d); continue
    if g not in cr.columns:
        d.update(status="无法检验", reason="不在 CRISPR 矩阵"); rows.append(d); continue
    x = cr[g].dropna(); xl = cr.loc[l_cr, g].dropna()
    frac_all = float((x <= -0.5).mean()); frac_l = float((xl <= -0.5).mean()) if len(xl) else np.nan
    growth = float((xl >= 0.3).mean()) if len(xl) else np.nan
    expr = float((ex.loc[l_ex, g] >= 1).mean()) if g in ex.columns and len(l_ex) else np.nan
    d.update(n_luad_crispr=len(xl), frac_dep_luad=round(frac_l, 3), frac_dep_all=round(frac_all, 3),
             median_chronos_luad=round(float(xl.median()), 3), median_chronos_nonluad=round(float(cr.loc[~cr.index.isin(luad_ids), g].median()), 3),
             frac_growth_luad=round(growth, 3), frac_expr_luad=round(expr, 3) if expr == expr else np.nan)
    reasons = []
    common = frac_all >= 0.9
    if common: reasons.append("常见必需基因（否决，长期毒性风险）")
    if growth >= 0.2: reasons.append("敲除促进生长（否决）")
    expressed = expr >= 0.5 if expr == expr else False
    dep = frac_l >= 0.2
    if not expressed: reasons.append("肺腺癌系中表达不足" if expr == expr else "无表达数据")
    if not dep: reasons.append("依赖比例<20%（无存活依赖）")
    ok = dep and expressed and not common and not growth >= 0.2
    d["common_essential"] = common
    d["status"] = "通过" if ok else "未通过"
    d["reason"] = "；".join(reasons)
    d["to_supplement_pool"] = (not ok) and (not common) and not (growth >= 0.2)
    rows.append(d)
out = pd.DataFrame(rows).set_index("gene")
out.to_csv(f"{P}/results/08_layer3_depmap/layer3_depmap.tsv", sep="\t")
up = out[out.direction == "up"]
print("上调基因", len(up), "；状态：", up["status"].value_counts().to_dict())
print("常见必需：", int(up["common_essential"].fillna(False).sum()), "；待补充池：", int(up["to_supplement_pool"].fillna(False).sum()))
for lab, sub in [("预设 352", up[up.in_main352]), ("683", up)]:
    print(lab, "上调", len(sub), sub["status"].value_counts().to_dict())
print("通过者（预设 352 内）：", ", ".join(f"{g}({r.frac_dep_luad:.2f})" for g, r in up[(up.status == '通过') & up.in_main352].iterrows()))
print("通过者（仅 AT2 谱系版）：", ", ".join(up[(up.status == '通过') & ~up.in_main352].index))
