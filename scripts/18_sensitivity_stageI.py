#!/usr/bin/env python
"""敏感性分析（Table S1 分期；方案在运行前写入 docs/项目进展与交接.md）：排除病理 II–III 期患者 P5、P11、P20、P24、P25（P7 仍排除），剩 17 对。
a) 上皮 A；b) 髓系 SPP1、TREM2、CXCL9；c) 增殖上皮比例配对比较。只报告，不改变任何入选判定。
用法：python scripts/18_sensitivity_stageI.py → results/18_sensitivity/"""
import os, sys, importlib.util, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import pearsonr, spearmanr, wilcoxon
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/18_sensitivity"; os.makedirs(OUT, exist_ok=True)
spec = importlib.util.spec_from_file_location("m14", f"{P}/scripts/14_microenv_layer1.py"); m14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m14)
EXCL = ["P5", "P11", "P20", "P24", "P25"]
def log(m): print(m, flush=True); open(f"{OUT}/run_log.txt", "a").write(m + "\n")

# ---------- a) 上皮 A ----------
cnt = pd.read_csv(f"{P}/results/03_pseudobulk/pseudobulk_counts_A_main.tsv", sep="\t", index_col=0)
meta = pd.read_csv(f"{P}/results/03_pseudobulk/pseudobulk_meta_A_main.tsv", sep="\t", index_col=0)
keep = ~meta["patient"].isin(EXCL)
cnt, meta = cnt[keep.values], meta[keep.values].copy()
meta["condition"] = pd.Categorical(meta["condition"], categories=["Precursor", "LUAD"])
log(f"[a] 上皮 A 敏感性：{len(meta)//2} 对；患者 {sorted(meta['patient'].unique(), key=lambda x:int(x[1:]))}")
r = m14.deseq(cnt, meta[["patient", "condition"]], "~ patient + condition", ["condition", "LUAD", "Precursor"])
r.to_csv(f"{OUT}/deg_A_stageI17.tsv", sep="\t")
main = pd.read_csv(f"{P}/results/03_pseudobulk/deg_A_main.tsv", sep="\t", index_col=0)
sig_s = r[r["padj"] < 0.05]; sig_m = main[main["padj"] < 0.05]
log(f"[a] 敏感性 检测基因 {int(r['padj'].notna().sum())}；显著 {len(sig_s)}（上调 {int((sig_s.log2FoldChange>0).sum())} / 下调 {int((sig_s.log2FoldChange<0).sum())}）；主分析显著 {len(sig_m)}")
common = r.index.intersection(main.index); cm = common[r.loc[common, "log2FoldChange"].notna() & main.loc[common, "log2FoldChange"].notna()]
x, y = main.loc[cm, "log2FoldChange"], r.loc[cm, "log2FoldChange"]
log(f"[a] 两者都检测到的基因 {len(cm)}；log2FC 相关 Pearson {pearsonr(x, y)[0]:.3f}，Spearman {spearmanr(x, y)[0]:.3f}")
ms = sig_m.index.intersection(r.index)
same = (np.sign(main.loc[ms, "log2FoldChange"]) == np.sign(r.loc[ms, "log2FoldChange"]))
log(f"[a] 主分析显著基因中两者都检测到 {len(ms)}；方向一致 {int(same.sum())}（{same.mean()*100:.1f}%）；其中敏感性也 FDR<0.05 且同向 {int((same & (r.loc[ms,'padj']<0.05)).sum())}")
c = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
c352 = c[c.in_main352]
t = r.reindex(c352.index)
ok = (t["padj"] < 0.05) & (np.sign(t["log2FoldChange"]) == np.sign(c352["log2FoldChange"]))
log(f"[a] 预设 352：敏感性中仍 FDR<0.05 且同向 {int(ok.sum())}/352（上调 {int(ok[c352.direction=='up'].sum())}/224，下调 {int(ok[c352.direction=='down'].sum())}/128）；有检测值 {int(t['padj'].notna().sum())}；方向一致 {int((np.sign(t['log2FoldChange'])==np.sign(c352['log2FoldChange'])).sum())}")
pd.DataFrame({"main_lfc": c352["log2FoldChange"], "sens_lfc": t["log2FoldChange"], "sens_padj": t["padj"], "retained": ok}).to_csv(f"{OUT}/layer1_352_retention.tsv", sep="\t")
# ---------- c) 增殖上皮比例 ----------
pt = pd.read_csv(f"{P}/results/03_pseudobulk/pairing_table.tsv", sep="\t")
def cmp(df, lab):
    d = (df["Prolif_frac_luad"] - df["Prolif_frac_precursor"]).values
    log(f"[c] {lab}：{len(df)} 对；前体均值 {df['Prolif_frac_precursor'].mean()*100:.1f}%，LUAD 均值 {df['Prolif_frac_luad'].mean()*100:.1f}%；均值差 {d.mean()*100:+.1f} 个百分点；中位差 {np.median(d)*100:+.1f}；LUAD 更高 {int((d>0).sum())}/{len(d)}；配对 Wilcoxon 双侧 P={wilcoxon(d).pvalue:.4f}")
    return d
cmp(pt, "主分析"); cmp(pt[~pt["patient"].isin(EXCL)], "敏感性（17 对）")
# 同样看 AT2/Club-like 与 AT2 构成
for col in ["AT2_frac", "ATClub_frac"]:
    for lab, df in [("主分析", pt), ("17 对", pt[~pt["patient"].isin(EXCL)])]:
        d = (df[f"{col}_luad"] - df[f"{col}_precursor"]).values
        log(f"[c] {col} {lab}：均值差 {d.mean()*100:+.1f} 个百分点；LUAD 更高 {int((d>0).sum())}/{len(d)}；Wilcoxon P={wilcoxon(d).pvalue:.4f}")
# ---------- b) 髓系 ----------
a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
a.obs["lin"] = ""
for L, labs in m14.LIN.items(): a.obs.loc[a.obs["final_label"].astype(str).isin(labs), "lin"] = L
epcam = np.asarray(a[:, "EPCAM"].layers["counts"].todense()).ravel()
a = a[((a.obs["lin"] == "Myeloid") & (epcam < 2)).values].copy()
a.obs["grp"] = a.obs["stage"].astype(str).map(lambda s: "LUAD" if s == "LUAD" else ("Precursor" if s in ("AAH", "AIS", "MIA") else "Normal"))
def myeloid(excl):
    groups, mt = {}, []
    for p in sorted(a.obs["patient"].unique(), key=lambda x: int(x[1:])):
        if p in excl: continue
        idxs = {g: np.where(((a.obs["patient"] == p) & (a.obs["grp"] == g)).values)[0] for g in ["Precursor", "LUAD"]}
        if min(len(v) for v in idxs.values()) < m14.MIN_CELLS: continue
        for g, v in idxs.items(): groups[f"{p}_{g}"] = v; mt.append(dict(sample=f"{p}_{g}", patient=p, condition=g))
    counts = pd.DataFrame(m14.sums(a, groups), index=a.var_names).T.astype(int)
    mt = pd.DataFrame(mt).set_index("sample"); mt["condition"] = pd.Categorical(mt["condition"], categories=["Precursor", "LUAD"])
    return m14.deseq(counts, mt[["patient", "condition"]], "~ patient + condition", ["condition", "LUAD", "Precursor"]), len(mt) // 2
rm, nm = myeloid([]); rs, ns = myeloid(EXCL)
rs.to_csv(f"{OUT}/deg_Myeloid_stageI.tsv", sep="\t")
log(f"[b] 髓系：主分析重算 {nm} 对；敏感性 {ns} 对")
for g in ["SPP1", "TREM2", "CXCL9"]:
    log(f"[b] {g}：主分析 log2FC={rm.loc[g,'log2FoldChange']:+.3f}，FDR={rm.loc[g,'padj']:.2e}；敏感性（{ns} 对）log2FC={rs.loc[g,'log2FoldChange']:+.3f}，FDR={rs.loc[g,'padj']:.2e}")
log("===== 完成 =====")
