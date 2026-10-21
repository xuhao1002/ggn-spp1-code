"""重构脚本（2026-10-21）：缺氧两基因集的细胞来源构成（稿件结果"post hoc check"一句的依据）。
原核查为交互式运行，脚本未保存；本脚本按交接文档记录的口径重构：每个基因取发现队列中表达中位数最高的谱系
（results/04_source_filter/ratio_peng.tsv 的 epi 与 med_* 列，与 98_niche_sensitivity.py / A6_prepare.py 的 best 定义相同），
统计 HIF13（HIF14_used）与 HALLMARK197（HALLMARK_used）中各谱系基因的比例。断言复现日志记录值：上皮 30% vs 26%，髓系 20% vs 15%。"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import pandas as pd
sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_hypoxia_composition")
G = json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
rat = pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv", sep="\t", index_col=0)
best = rat[["epi"] + [c for c in rat.columns if c.startswith("med_")]].idxmax(axis=1)
res = {"plan_sha256": sha, "source": "results/04_source_filter/ratio_peng.tsv (top lineage by median expression)"}
for nm, lab in [("HIF14_used", "HIF13"), ("HALLMARK_used", "HALLMARK197")]:
    gs = [g for g in G[nm] if g in best.index]
    vc = best[gs].value_counts(normalize=True)
    res[lab] = {"n_genes_in_set": len(G[nm]), "n_genes_with_lineage": len(gs), "fractions": {k: float(v) for k, v in vc.items()},
                "missing_from_table": sorted(set(G[nm]) - set(gs))}
    pd.DataFrame({"gene": gs, "top_lineage": best[gs].values}).to_csv(OUT + f"composition_{lab}.tsv", sep="\t", index=False)
h, a = res["HIF13"]["fractions"], res["HALLMARK197"]["fractions"]
assert round(h["epi"], 2) == 0.30 and round(a["epi"], 2) == 0.26 and round(h["med_Myeloid"], 2) == 0.20 and round(a["med_Myeloid"], 2) == 0.15, (h, a)
res["reproduces_logged_values"] = True
dump(res, OUT + "composition_results.json")
log("完成", {k: round(v, 3) for k, v in h.items()}, {k: round(v, 3) for k, v in a.items()})
