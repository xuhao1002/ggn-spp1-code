#!/usr/bin/env python
"""间接路径（标准见 docs/项目进展与交接.md，2026-10-08 用户确认）。
- TF 候选：预设 352 内（AT2 谱系版并列）上调基因，DGIdb 无抑制类药物（第 5 层 drug_grade=='无'），且在调控资源中有 >=5 个靶基因。
- 调控资源：CollecTRI（decoupler 2.2.0，2026-10-08 取自 OmniPath，42,990 条）∪ DoRothEA A/B 级（15,117 条）；两者符号冲突的关系丢弃。
- 效应分子：必须完整满足第 1 层（预设 352 版；AT2 谱系版并列），方向=TF 方向×调控符号，通过来源过滤，并独立通过第 3–5 层（第 3 层=DepMap 通过 或 LINCS 通过）。
- TF 活性（ULM，伪批量 A 的 22 对，配对检验，BH）只作辅助。
输出 results/13_indirect/"""
import os, numpy as np, pandas as pd, decoupler as dc
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/13_indirect"
I = pd.read_csv(f"{P}/results/11_integrated/integrated_layers.tsv", sep="\t", index_col=0)
l3l = pd.read_csv(f"{P}/results/12_layer3_lincs/layer3_lincs.tsv", sep="\t", index_col=0)
I["L3_lincs"] = l3l["status"].reindex(I.index)
I["L3_any"] = (I["L3"] == "通过") | (I["L3_lincs"] == "通过")
col = pd.read_csv(f"{P}/raw/annotation/collectri_human.tsv", sep="\t"); dor = pd.read_csv(f"{P}/raw/annotation/dorothea_AB_human.tsv", sep="\t")
col["sign"] = np.sign(col["weight"]); dor["sign"] = np.sign(dor["weight"])
net = pd.concat([col[["source", "target", "sign"]].assign(res="CollecTRI"), dor[["source", "target", "sign"]].assign(res="DoRothEA_AB")])
g = net.groupby(["source", "target"])["sign"].agg(["min", "max"]).reset_index()
bad = g[g["min"] != g["max"]]; print("符号冲突关系", len(bad))
net = net.merge(g[g["min"] == g["max"]][["source", "target"]], on=["source", "target"]).drop_duplicates(["source", "target"])
print("合并后关系", len(net), "；TF", net["source"].nunique())
ntar = net.groupby("source").size()

def run(scope_col, label):
    sub = I[(I.direction == "up") & I[scope_col]]
    tf = sub[(sub["L5_drug_grade"] == "无") & sub.index.isin(ntar[ntar >= 5].index)]
    print(f"[{label}] 上调 {len(sub)}；无抑制药物 {int((sub['L5_drug_grade']=='无').sum())}；其中为 TF（>=5 靶基因）{len(tf)}")
    rows = []
    for f, r in tf.iterrows():
        edges = net[net["source"] == f]
        for _, e in edges.iterrows():
            t = e["target"]
            if t not in I.index: continue
            x = I.loc[t]
            exp_dir = np.sign(r["log2FoldChange"]) * e["sign"]
            ok_dir = np.sign(x["log2FoldChange"]) == exp_dir
            l1 = bool(x[scope_col]) and ok_dir                     # 第 1 层完整满足（该版候选池内、方向一致）
            rows.append(dict(TF=f, TF_lfc=round(r["log2FoldChange"], 2), target=t, sign=int(e["sign"]), res=e["res"], target_lfc=round(x["log2FoldChange"], 2),
                             L1=l1, L3=bool(x["L3_any"]), L4=x["L4"] == "通过", L5=x["L5"] == "通过", target_direction=x["direction"],
                             drug_grade=x["L5_drug_grade"], approved=x["L5_approved_drugs"]))
    d = pd.DataFrame(rows)
    return tf, d
res = {}
for scope, label in [("in_main352", "预设352"), ("in_any", "AT2谱系版并列（683）")]:
    if scope == "in_any": I["in_any"] = True
    tf, d = run(scope, label)
    d["all_pass"] = d[["L1", "L3", "L4", "L5"]].all(axis=1) if len(d) else []
    if len(d):
        eff = d[d.all_pass & (d.target_direction == "up") | d.all_pass]
        print(f"[{label}] TF {len(tf)}；TF→靶基因关系（靶基因在候选表内）{len(d)}；满足第1层 {int(d.L1.sum())}；独立通过第3–5层 {int((d.L1 & d.L3 & d.L4 & d.L5).sum())}；效应分子（全部满足）{int(d.all_pass.sum())}")
        print("  按层通过（在满足第1层的关系中）：L3", int((d.L1 & d.L3).sum()), "L4", int((d.L1 & d.L4).sum()), "L5", int((d.L1 & d.L5).sum()),
              "；L4&L5", int((d.L1 & d.L4 & d.L5).sum()), "；L3&L4", int((d.L1 & d.L3 & d.L4).sum()))
        print("  效应分子：", d[d.all_pass][["TF", "target"]].drop_duplicates().to_dict("records"))
    d.to_csv(f"{OUT}/tf_effectors_{label.split('（')[0]}.tsv", sep="\t", index=False)
    tf.to_csv(f"{OUT}/tf_candidates_{label.split('（')[0]}.tsv", sep="\t")
    res[label] = (tf, d)
tf352 = res["预设352"][0]
print("预设352 内 TF 候选：", ", ".join(f"{g}(靶基因{int(ntar[g])})" for g in tf352.index))

# TF 活性（ULM，伪批量 A，22 对）——只报告
cnt = pd.read_csv(f"{P}/results/03_pseudobulk/pseudobulk_counts_A_main.tsv", sep="\t", index_col=0)
lib = cnt.div(cnt.sum(axis=1), axis=0) * 1e6; X = np.log1p(lib)
tfs = [t for t in tf352.index if t in ntar.index]
netu = net.rename(columns={"sign": "weight"})[["source", "target", "weight"]]
netu = netu[netu["source"].isin(tfs)]
est = dc.mt.ulm(data=X, net=netu, tmin=5)
act = est[0] if isinstance(est, tuple) else est
rows = []
for f in act.columns:
    pre = act.loc[[i for i in act.index if i.endswith("_Precursor")], f]; lu = act.loc[[i for i in act.index if i.endswith("_LUAD")], f]
    pre.index = [i.split("_")[0] for i in pre.index]; lu.index = [i.split("_")[0] for i in lu.index]
    dd = (lu - pre.reindex(lu.index)).dropna()
    if len(dd) >= 5 and dd.abs().sum() > 0:
        rows.append(dict(TF=f, n_pairs=len(dd), mean_diff=round(dd.mean(), 3), p=wilcoxon(dd).pvalue, n_up=int((dd > 0).sum())))
ta = pd.DataFrame(rows).set_index("TF")
if len(ta):
    ta["q"] = multipletests(ta["p"], method="fdr_bh")[1]; ta.to_csv(f"{OUT}/tf_activity_ULM.tsv", sep="\t")
    print("TF 活性（ULM，配对 Wilcoxon，BH）：", len(ta), "个 TF；q<0.05 且 LUAD 更高：", ", ".join(ta[(ta.q < 0.05) & (ta.mean_diff > 0)].index))
