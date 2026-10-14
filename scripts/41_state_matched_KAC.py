#!/usr/bin/env python
"""第四部分 第 2 步（缩减版，按方案 2.3 与 §10）：KAC 可配对仅 7 对 <8，KAC 比较「无法进行」，不下调门槛；
只做 AT2L-other（KAC 得分 <= 固定阈值 0.1196 的 AT2 谱系细胞）前体 vs Normal 配对伪批量（pydeseq2 ~ patient + condition，每侧>=50 细胞）。
用途：检验第二阶段 AT2L 候选的升高是否在非 KAC 细胞中仍存在（即不只是 KAC 细胞比例变化所致）。不分梯队 S/R（缺 KAC 侧比较），不产生新候选。
本研究不推断干预方向的正负。输出 results/40_drivers/state_matched/"""
import os, re, sys, time, importlib.util, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/40_drivers/state_matched"
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("m14", f"{P}/scripts/14_microenv_layer1.py"); m14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m14)
LOG = open(f"{P}/results/40_drivers/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()
txt = open(f"{P}/results/20_stage2/discovery/run_log.txt", encoding="utf-8").read()
THR = float(re.search(r"\{'KAC': np\.float64\(([0-9.]+)\)", txt).group(1))
AT2L = ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"]
kac = [l.strip() for l in open(f"{P}/results/20_stage2/peng_TableS3_KAC_human.txt")]
log(f"===== 第 2 步（缩减）：AT2L-other 对照；KAC 比较无法进行（可配对 7 对 <8）；阈值 {THR:.6f} =====")
a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad"); o = a.obs
o["patient"] = o["patient"].astype(str); o["stage"] = o["stage"].astype(str)
m = o["final_label"].astype(str).isin(AT2L).values
s2 = a[m].copy(); sc.tl.score_genes(s2, [g for g in kac if g in s2.var_names], score_name="KAC", random_state=0)
rec = np.percentile(s2.obs.loc[s2.obs.stage == "Normal", "KAC"], 90); log(f"重现性核对：差 {abs(rec-THR):.2e}")
s2.obs["other"] = s2.obs["KAC"].values <= THR
X = s2.layers["counts"].tocsr(); ob = s2.obs; idxs, meta = {}, []
for p in sorted(ob.patient.unique(), key=lambda x: int(x[1:])):
    for g, st in [("Normal", ["Normal"]), ("Precursor", ["AAH", "AIS", "MIA"])]:
        ix = np.where((ob.patient.values == p) & np.isin(ob.stage.values, st) & ob.other.values)[0]; idxs[f"{p}_{g}"] = ix
        meta.append(dict(sample=f"{p}_{g}", patient=p, condition=g, n_cells=len(ix)))
meta = pd.DataFrame(meta).set_index("sample"); ok = meta.groupby("patient").n_cells.min() >= 50; keep = ok[ok].index
meta = meta[meta.patient.isin(keep)]; log(f"AT2L-other 可配对 {len(keep)} 对；排除：{list(ok[~ok].index)}")
counts = pd.DataFrame({k: np.asarray(X[idxs[k]].sum(axis=0)).ravel() for k in meta.index}, index=s2.var_names).T.astype(int)
mt = meta.copy(); mt["condition"] = pd.Categorical(mt["condition"], categories=["Normal", "Precursor"])
r = m14.deseq(counts, mt[["patient", "condition"]], "~ patient + condition", ["condition", "Precursor", "Normal"])
r["p1_up"] = np.where(r["log2FoldChange"] > 0, r["pvalue"] / 2, 1 - r["pvalue"] / 2); r.to_csv(f"{OUT}/deg_AT2L_other.tsv", sep="\t"); meta.to_csv(f"{OUT}/units_AT2L_other.tsv", sep="\t")
log(f"AT2L-other：检测基因 {int(r.padj.notna().sum())}；FDR<0.05 上调 {int(((r.padj<0.05)&(r.log2FoldChange>0)).sum())}，下调 {int(((r.padj<0.05)&(r.log2FoldChange<0)).sum())}")
fin = pd.read_csv(f"{P}/results/20_stage2/stage2_pass_all.tsv", sep="\t")
main = pd.read_csv(f"{P}/results/20_stage2/discovery/deg_AT2L_main.tsv", sep="\t", index_col=0)
fin["main_lfc"] = [main.log2FoldChange.get(g, np.nan) if l == "AT2L" else np.nan for g, l in zip(fin.gene, fin.lineage)]
fin["other_lfc"] = [r.log2FoldChange.get(g, np.nan) if l == "AT2L" else np.nan for g, l in zip(fin.gene, fin.lineage)]
fin["other_padj"] = [r.padj.get(g, np.nan) if l == "AT2L" else np.nan for g, l in zip(fin.gene, fin.lineage)]
fin["persists_in_other"] = (fin.other_lfc > 0) & (fin.other_padj < 0.05)
fin["KAC_comparison"] = "无法进行（可配对 7 对 <8）"; fin["note"] = "本研究不推断干预方向的正负"
cols = [c for c in ["gene", "lineage", "main_lfc", "other_lfc", "other_padj", "persists_in_other", "KAC_comparison", "note"] if c in fin.columns]
fin[cols].to_csv(f"{OUT}/stage2_32genes_state_matched.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); log("32 个基因×群在 AT2L-other 中：\n" + fin[cols[:6]].round(4).to_string(index=False))
log(f"AT2L 的 31 个基因中，在 AT2L-other 仍上调且 FDR<0.05：{int(fin.persists_in_other.sum())}；方向仍为正：{int((fin.other_lfc>0).sum())}")
log("===== 第 2 步（缩减）完成 =====")
