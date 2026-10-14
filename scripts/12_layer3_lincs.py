#!/usr/bin/env python
"""第 3 层补充：LINCS L1000 CRISPR KO 共识特征（标准见 docs/项目进展与交接.md「补充分析判定标准」，2026-10-08 用户确认）。
- 程序：预设 352 内经验证的上调 224 + 下调 128，只保留 L1000 测到的基因（gene_info 中 pr_is_lm==1 或 pr_is_bing==1）；
  拆为增殖（Tirosh 96 ∪ GO 细胞周期相关）与非增殖。
- 每个池内基因的 KO 特征（上调集 U、下调集 D，限于测到的基因）与程序的重叠用超几何检验（背景=测到的基因数）。
  S = [−log10 p(D∩程序上调) − log10 p(U∩程序下调)] − [−log10 p(U∩程序上调) − log10 p(D∩程序下调)]（正值=敲除后程序被逆转）。
- 零分布：库内其他基因 KO 特征的 S（程序基因自身作为 KO 对象的不入零分布）；经验单侧 P=(1+#null>=S)/(1+N)；池内基因 BH 校正（非增殖、增殖分别校正）。
- 通过=非增殖程序 q<0.05 且 S>0；仅增殖部分显著逆转=「仅增殖逆转」，不算通过。CRISPR 不要求目标基因自身下降。
输出 results/12_layer3_lincs/layer3_lincs.tsv"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.expanduser("~/全部数据库/GGN_target/scripts"))
import prolif_genes
P = os.path.expanduser("~/全部数据库/GGN_target"); R = f"{P}/raw/lincs"; OUT = f"{P}/results/12_layer3_lincs"
gi = pd.read_csv(f"{R}/gene_info.txt.gz", sep="\t", compression="gzip", dtype=str)
measured = set(gi.loc[(gi["pr_is_lm"] == "1") | (gi["pr_is_bing"] == "1"), "pr_gene_symbol"])
M = len(measured)
cand = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
l3 = pd.read_csv(f"{P}/results/08_layer3_depmap/layer3_depmap.tsv", sep="\t", index_col=0)
pool = l3[(l3.direction == "up") & (l3["to_supplement_pool"] == True)].index.tolist()
prolif = prolif_genes.prolif_set()
main = cand[cand.in_main352]
up_all, dn_all = set(main.index[main.direction == "up"]), set(main.index[main.direction == "down"])
up_m, dn_m = up_all & measured, dn_all & measured
prog = {"np_up": up_m - prolif, "np_dn": dn_m - prolif, "pr_up": up_m & prolif, "pr_dn": dn_m & prolif}
print(f"L1000 测到的基因 {M}；程序上调 {len(up_all)}→{len(up_m)}、下调 {len(dn_all)}→{len(dn_m)}；非增殖 上 {len(prog['np_up'])} 下 {len(prog['np_dn'])}；增殖 上 {len(prog['pr_up'])} 下 {len(prog['pr_dn'])}")

# 读 KO 库
sig = {}
for line in open(f"{R}/Enrichr_LINCS_L1000_CRISPR_KO_Consensus_Sigs.txt"):
    parts = line.rstrip("\n").split("\t"); name = parts[0]; genes = {g for g in parts[2:] if g}
    if name.endswith(" Up"): sig.setdefault(name[:-3], {})["U"] = genes & measured
    elif name.endswith(" Down"): sig.setdefault(name[:-5], {})["D"] = genes & measured
sig = {k: v for k, v in sig.items() if "U" in v and "D" in v}
print("KO 基因数", len(sig))

def lp(k, n_prog, n_sig):
    if n_prog == 0 or n_sig == 0: return 0.0
    return float(-np.log10(max(hypergeom.sf(k - 1, M, n_prog, n_sig), 1e-300)))

def S_of(U, D, up, dn):
    rev = lp(len(D & up), len(up), len(D)) + lp(len(U & dn), len(dn), len(U))
    con = lp(len(U & up), len(up), len(U)) + lp(len(D & dn), len(dn), len(D))
    return rev - con

prog_genes = up_all | dn_all
null_genes = [g for g in sig if g not in prog_genes]
null_np = np.array([S_of(sig[g]["U"], sig[g]["D"], prog["np_up"], prog["np_dn"]) for g in null_genes])
null_pr = np.array([S_of(sig[g]["U"], sig[g]["D"], prog["pr_up"], prog["pr_dn"]) for g in null_genes])
print("零分布基因数", len(null_genes))
rows = []
for g in pool:
    d = dict(gene=g, set=cand.loc[g, "set"], in_main352=cand.loc[g, "in_main352"], tier=cand.loc[g, "tier"])
    if g not in sig:
        d["status"] = "无法检验"; d["reason"] = "不在 KO 共识库"; rows.append(d); continue
    U, D = sig[g]["U"], sig[g]["D"]
    s_np = S_of(U, D, prog["np_up"], prog["np_dn"]); s_pr = S_of(U, D, prog["pr_up"], prog["pr_dn"])
    d.update(S_nonprolif=round(s_np, 3), S_prolif=round(s_pr, 3),
             p_nonprolif=(1 + (null_np >= s_np).sum()) / (1 + len(null_np)), p_prolif=(1 + (null_pr >= s_pr).sum()) / (1 + len(null_pr)),
             n_U=len(U), n_D=len(D), self_in_D=g in D)
    rows.append(d)
out = pd.DataFrame(rows).set_index("gene")
t = out["p_nonprolif"].notna()
out.loc[t, "q_nonprolif"] = multipletests(out.loc[t, "p_nonprolif"], method="fdr_bh")[1]
out.loc[t, "q_prolif"] = multipletests(out.loc[t, "p_prolif"], method="fdr_bh")[1]
def st(r):
    if r["status"] == "无法检验" if "status" in r and isinstance(r["status"], str) else False: return "无法检验"
    if r["q_nonprolif"] < 0.05 and r["S_nonprolif"] > 0: return "通过"
    if r["q_prolif"] < 0.05 and r["S_prolif"] > 0: return "仅增殖逆转"
    return "未通过"
out["status"] = [("无法检验" if not t[g] else st(out.loc[g])) for g in out.index]
out.to_csv(f"{OUT}/layer3_lincs.tsv", sep="\t")
print("池内基因", len(out), "；覆盖（在 KO 库中）", int(t.sum()), "（预设 352 内 %d/%d）" % (int((t & out.in_main352).sum()), int(out.in_main352.sum())))
print(out["status"].value_counts().to_dict(), "；预设 352 内：", out[out.in_main352]["status"].value_counts().to_dict())
print("通过者：", ", ".join(f"{g}(S={r.S_nonprolif:.1f},q={r.q_nonprolif:.3f})" for g, r in out[out.status == '通过'].sort_values('q_nonprolif').iterrows()))
print("仅增殖逆转：", ", ".join(out[out.status == '仅增殖逆转'].index))
