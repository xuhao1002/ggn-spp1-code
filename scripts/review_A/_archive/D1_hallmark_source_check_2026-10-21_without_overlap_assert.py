"""偏离记录 D1 结案的可复现脚本：比较发现队列所用 HALLMARK 低氧基因集（Enrichr MSigDB_Hallmark_2020 "Hypoxia"）
与 MSigDB v2023.2.Hs 官方 HALLMARK_HYPOXIA。官方 gmt 于 2026-10-21 自
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt 下载并存档于
results/review_D1/run_20261021_1400/（首次比较在云端容器完成）；本脚本离线读取存档文件并断言其 SHA-256。
所用列表 = results/hypoxia/gene_sets_used.json 的 HALLMARK_used（197）∪ 与 30 基因集合重叠而被移除的基因。"""
import os, sys, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import pandas as pd
sha = check_plan(); log("方案哈希校验通过", sha[:16])
GMT = f"{P}/results/review_D1/run_20261021_1400/h.all.v2023.2.Hs.symbols.gmt"
GMT_SHA = "8ff1f036f0988c99d42cf9e5ea24d355f10b8ece5e0d4f6edf26a0458f480780"
assert hashlib.sha256(open(GMT, "rb").read()).hexdigest() == GMT_SHA, "官方 gmt 存档的 SHA-256 不符"
OUT = outdir("review_D1")
off = [l.rstrip("\n").split("\t")[2:] for l in open(GMT) if l.startswith("HALLMARK_HYPOXIA\t")][0]
G = json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
top = pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv", sep="\t").gene.tolist()
removed = sorted(set(top) & set(off))
used = set(G["HALLMARK_used"]) | set(removed)
assert len(G["HALLMARK_used"]) == 197 and set(G["HALLMARK_used"]).isdisjoint(top)
res = {"plan_sha256": sha, "official_gmt_sha256": GMT_SHA, "n_official": len(off), "n_used_reconstructed": len(used),
       "removed_overlap_with_30_gene_set": removed, "only_used": sorted(used - set(off)), "only_official": sorted(set(off) - used),
       "identical": used == set(off),
       "note": "Reconstruction assumes the Enrichr list equals HALLMARK_used plus its overlap with the 30-gene set; A6c asserted the Enrichr list has 200 genes and that HALLMARK_used equals it minus that overlap."}
dump(res, OUT + "D1_HALLMARK_comparison.json")
log("D1：完全相同" if res["identical"] else "D1：不同", res["only_used"], res["only_official"])
