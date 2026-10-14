#!/usr/bin/env python
"""第四部分 第 4 步（分析 D）：H4 TIM-3/HAVCR2 假说（预设见方案第 4.1 节，用户已确认）。只报告，不入选。
- H4a T/NK HAVCR2 伪批量表达；H4b 髓系 HAVCR2；H4c T/NK HAVCR2+ 细胞比例（阈值 = Normal T/NK 细胞 HAVCR2 log1p 归一化表达第 90 百分位，本次运行一次性计算并固定；
  若该百分位为 0，则按第二阶段 IL1B 的同一规则取 HAVCR2>0，并如实记录）；H4d AT2 谱系 LGALS9、CEACAM1 表达（各自检验）。
- 分步（AAH、AIS/MIA）各报告；每步 H4 家族共 5 个检验，Holm（家族=5 个检验），同时报告未校正 P。表达检验：pydeseq2 `~ patient + condition` 单侧；单细胞最低 8 对。
- 验证：GSE164789（AIS/MIA 病灶 vs 同患者癌旁，可配对 <8 对标「无法验证」，阈值取 -A 细胞第 90 百分位）；GSE288479 bulk（7 对 GG vs DN；表达类，log2 CPM 配对 Wilcoxon 单侧）。
- 文献：Zhu B 等 Cancer Cell 2025;43:1125（PMID 40345189），原文称 putative target，不得写成已证实的人体 GGN 预防药物。早期阻断与较晚干预效果不同。
本研究不推断干预方向的正负。输出 results/40_drivers/hypotheses/"""
import os, sys, time, importlib.util, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/40_drivers/hypotheses"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("m14", f"{P}/scripts/14_microenv_layer1.py"); m14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m14)
LOG = open(f"{P}/results/40_drivers/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()
NOTE = "本研究不推断干预方向的正负"
MIN_CELLS, MIN_PAIRS = 50, 8
log("===== 第四部分 第 4 步 分析 D（H4：HAVCR2/TIM-3）=====")
def wil(d):
    d = np.asarray(d, float); d = d[~np.isnan(d)]
    if len(d) < 5 or np.all(d == 0): return np.nan
    return float(wilcoxon(d, alternative="greater").pvalue)
def p1_of(r): return np.where(r["log2FoldChange"] > 0, r["pvalue"] / 2, 1 - r["pvalue"] / 2)
def dense(M): return np.asarray(M.todense()).ravel() if sp.issparse(M) else np.asarray(M).ravel()

# ================= Peng snRNA =================
LINS = {"T/NK": ["T cell", "NK", "Proliferating T/NK"], "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"], "AT2L": ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"]}
a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
o = a.obs; o["patient"] = o["patient"].astype(str); o["stage"] = o["stage"].astype(str); o["fl"] = o["final_label"].astype(str)
epcam = dense(a[:, "EPCAM"].layers["counts"]); X = a.layers["counts"].tocsr(); hav = dense(a[:, "HAVCR2"].X)
patients = sorted(o["patient"].unique(), key=lambda x: int(x[1:]))
STEPS = {"AAH": ["AAH"], "AISMIA": ["AIS", "MIA"]}
# H4c 阈值：一次性计算并固定
m_t = o["fl"].isin(LINS["T/NK"]).values & (epcam < 2)
t_nrm = hav[m_t & (o["stage"].values == "Normal")]
p90 = float(np.percentile(t_nrm, 90)); rule = "> 第90百分位" if p90 > 0 else "HAVCR2>0（第 90 百分位为 0，按 IL1B 同规则）"
THR_H = p90
HPOS = (lambda v: v > THR_H) if p90 > 0 else (lambda v: v > 0)
log(f"H4c 阈值（一次性固定）：Normal T/NK 细胞（{len(t_nrm)} 个）HAVCR2 log1p 归一化表达第 90 百分位 = {p90:.6f}；定义 {rule}；Normal 中 HAVCR2+ 比例 {HPOS(t_nrm).mean()*100:.1f}%")
open(f"{OUT}/H4c_threshold_fixed.txt", "w").write(f"HAVCR2 P90 (Normal T/NK) = {p90:.10f}\nrule = {rule}\n")

def units(labs, epf, step):
    m_ = o["fl"].isin(labs).values & ((epcam < 2) if epf else True); idxs, meta = {}, []
    for p in patients:
        sp_ = o["patient"].values == p
        nrm = np.where(m_ & sp_ & (o["stage"].values == "Normal"))[0]; pre = np.where(m_ & sp_ & np.isin(o["stage"].values, STEPS[step]))[0]
        if len(nrm) < MIN_CELLS or len(pre) < MIN_CELLS: continue
        for g, ix in [("Normal", nrm), ("Precursor", pre)]: idxs[f"{p}_{g}"] = ix; meta.append(dict(sample=f"{p}_{g}", patient=p, condition=g, n_cells=len(ix)))
    return idxs, (pd.DataFrame(meta).set_index("sample") if meta else pd.DataFrame())
tests = []
for step in STEPS:
    for tid, lin, gene in [("H4a", "T/NK", "HAVCR2"), ("H4b", "Myeloid", "HAVCR2"), ("H4d-LGALS9", "AT2L", "LGALS9"), ("H4d-CEACAM1", "AT2L", "CEACAM1")]:
        idxs, meta = units(LINS[lin], lin != "AT2L", step); npairs = len(meta) // 2
        if npairs < MIN_PAIRS: tests.append(dict(cohort="Peng", step=step, test=tid, lineage=lin, gene=gene, n_pairs=npairs, status=f"无法进行（{npairs} 对）")); continue
        counts = pd.DataFrame({k: np.asarray(X[ix].sum(axis=0)).ravel() for k, ix in idxs.items()}, index=a.var_names).T.astype(int)
        mt = meta.copy(); mt["condition"] = pd.Categorical(mt["condition"], categories=["Normal", "Precursor"])
        r = m14.deseq(counts, mt[["patient", "condition"]], "~ patient + condition", ["condition", "Precursor", "Normal"]); r["p1_up"] = p1_of(r)
        if gene not in r.index or pd.isna(r.loc[gene, "pvalue"]): tests.append(dict(cohort="Peng", step=step, test=tid, lineage=lin, gene=gene, n_pairs=npairs, status="过滤后无检测值")); continue
        lg = np.log1p(counts.div(counts.sum(axis=1), axis=0) * 1e6); pts = list(meta.patient.unique())
        dif = pd.Series({p: lg.loc[f"{p}_Precursor", gene] - lg.loc[f"{p}_Normal", gene] for p in pts})
        tests.append(dict(cohort="Peng", step=step, test=tid, lineage=lin, gene=gene, n_pairs=npairs, statistic=float(r.loc[gene, "log2FoldChange"]), p_one_sided=float(r.loc[gene, "p1_up"]),
                          p_two_sided=float(r.loc[gene, "pvalue"]), wilcoxon_p1=wil(dif), n_up=int((dif > 0).sum()), status="已检验", metric="pseudobulk log2FC"))
    idxs, meta = units(LINS["T/NK"], True, step); npairs = len(meta) // 2
    if npairs < MIN_PAIRS: tests.append(dict(cohort="Peng", step=step, test="H4c", lineage="T/NK", gene="HAVCR2+ 比例", n_pairs=npairs, status=f"无法进行（{npairs} 对）"))
    else:
        pts = list(meta.patient.unique())
        fr = pd.Series({k: float(HPOS(hav[ix]).mean()) for k, ix in idxs.items()}); dd = pd.Series({p: fr[f"{p}_Precursor"] - fr[f"{p}_Normal"] for p in pts})
        tests.append(dict(cohort="Peng", step=step, test="H4c", lineage="T/NK", gene="HAVCR2+ 比例", n_pairs=npairs, statistic=float(dd.mean()), p_one_sided=wil(dd), n_up=int((dd > 0).sum()), status="已检验", metric="mean fraction diff"))
    log(f"[{step}] 完成")
del a, X
T = pd.DataFrame(tests); T["p_holm_H4_step"] = np.nan
for step, g in T.groupby("step"):
    ok = g["p_one_sided"].notna()
    if ok.sum(): T.loc[g.index[ok.values], "p_holm_H4_step"] = multipletests(g.loc[ok, "p_one_sided"], method="holm")[1]
T["supported"] = (T["statistic"] > 0) & (T["p_one_sided"] < 0.05) & (T["p_holm_H4_step"] < 0.05)
T["note"] = NOTE + "；TIM-3 为 putative target"; T.to_csv(f"{OUT}/H4_peng_stepwise.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
log("H4 Peng 分步（" + NOTE + "）：\n" + T[["step", "test", "lineage", "gene", "n_pairs", "statistic", "p_one_sided", "p_holm_H4_step", "supported", "wilcoxon_p1", "n_up", "status"]].round(5).to_string(index=False))

# ================= GSE164789 =================
VL = {"T/NK": ["T cell", "NK"], "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"], "AT2L": ["AT2", "Club/secretory"]}
EPI_V = ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal"]
d = f"{P}/results/05_validation/GSE164789"
v = sc.read_h5ad(f"{d}/GSE164789_annotated.h5ad"); ann = pd.read_csv(f"{d}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
v.obs["fl"] = v.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
for c in ["patient", "tissue", "histology", "sample"]: v.obs[c] = v.obs[c].astype(str)
ep = dense(v[:, "EPCAM"].layers["counts"]); pt = dense(v[:, "PTPRC"].layers["counts"]); isepi = v.obs["fl"].isin(EPI_V).values
clean = ~(isepi & (pt >= 2)) & ~(~isepi & (ep >= 2)); v = v[clean & v.obs["histology"].isin(["AIS", "MIA"]).values].copy()
v.obs["unit"] = v.obs["patient"] + "_" + v.obs["histology"]; Xv = v.layers["counts"].tocsr(); hv = dense(v[:, "HAVCR2"].X)
def v_units(labs):
    m_ = v.obs["fl"].isin(labs).values; idxs, meta = {}, []
    for u in sorted(v.obs["unit"].unique()):
        a_ = np.where(m_ & (v.obs["unit"].values == u) & (v.obs["tissue"].values == "A"))[0]; t_ = np.where(m_ & (v.obs["unit"].values == u) & (v.obs["tissue"].values == "T"))[0]
        if len(a_) < 50 or len(t_) < 50: continue
        for g, ix in [("A", a_), ("T", t_)]: idxs[f"{u}|{g}"] = ix; meta.append(dict(sample=f"{u}|{g}", pair=u, tissue=g, n_cells=len(ix)))
    return idxs, (pd.DataFrame(meta).set_index("sample") if meta else pd.DataFrame())
vt = []
for tid, lin, gene in [("H4a", "T/NK", "HAVCR2"), ("H4b", "Myeloid", "HAVCR2"), ("H4d-LGALS9", "AT2L", "LGALS9"), ("H4d-CEACAM1", "AT2L", "CEACAM1"), ("H4c", "T/NK", "HAVCR2+ 比例")]:
    idxs, meta = v_units(VL[lin]); npairs = len(meta) // 2
    if npairs < MIN_PAIRS: vt.append(dict(cohort="GSE164789", test=tid, lineage=lin, gene=gene, n_pairs=npairs, status=f"无法验证（可配对 {npairs} 对 <{MIN_PAIRS}）")); continue
    if tid == "H4c":
        # 阈值取 -A 细胞第 90 百分位（与第二阶段做法一致）
        m_a = v.obs["fl"].isin(VL[lin]).values & (v.obs["tissue"].values == "A"); pa90 = np.percentile(hv[m_a], 90); pos = (lambda x: x > pa90) if pa90 > 0 else (lambda x: x > 0)
        fr = pd.Series({k: float(pos(hv[ix]).mean()) for k, ix in idxs.items()}); pairs = meta.pair.unique(); dd = pd.Series({u: fr[f"{u}|T"] - fr[f"{u}|A"] for u in pairs})
        vt.append(dict(cohort="GSE164789", test=tid, lineage=lin, gene=gene, n_pairs=npairs, statistic=float(dd.mean()), p_one_sided=wil(dd), n_up=int((dd > 0).sum()), status="已检验")); continue
    counts = pd.DataFrame({k: np.asarray(Xv[ix].sum(axis=0)).ravel() for k, ix in idxs.items()}, index=v.var_names).T.astype(int)
    mt = meta.copy(); mt["tissue"] = pd.Categorical(mt["tissue"], categories=["A", "T"])
    r = m14.deseq(counts, mt[["pair", "tissue"]], "~ pair + tissue", ["tissue", "T", "A"]); r["p1_up"] = p1_of(r)
    if gene in r.index and pd.notna(r.loc[gene, "pvalue"]): vt.append(dict(cohort="GSE164789", test=tid, lineage=lin, gene=gene, n_pairs=npairs, statistic=float(r.loc[gene, "log2FoldChange"]), p_one_sided=float(r.loc[gene, "p1_up"]), status="已检验"))
    else: vt.append(dict(cohort="GSE164789", test=tid, lineage=lin, gene=gene, n_pairs=npairs, status="过滤后无检测值"))
VT = pd.DataFrame(vt)
if VT["p_one_sided"].notna().any():
    ok = VT["p_one_sided"].notna(); VT.loc[ok, "p_holm_H4"] = multipletests(VT.loc[ok, "p_one_sided"], method="holm")[1]
else: VT["p_holm_H4"] = np.nan
log("H4 在 GSE164789（AIS/MIA 病灶 vs 癌旁）：\n" + VT.round(5).to_string(index=False))

# ================= GSE288479 bulk =================
BULK7 = ["SMC06", "SMC08", "SMC11", "SMC13", "SMC18", "SMC19", "SMC20"]
bk = pd.read_csv(f"{P}/raw/multi_GSE288479_partsolid/GSE288479_raw_counts_All_samples.txt.gz", index_col=0, compression="gzip")
bk = bk.drop(columns=["gene_id", "gene_type"]).groupby("gene_name").sum(); bk = bk[[f"{p}_{s}_counts" for p in BULK7 for s in ["DN", "GG"]]]
lb = np.log2(bk.div(bk.sum(axis=0), axis=1) * 1e6 + 1); bt = []
for g in ["HAVCR2", "LGALS9", "CEACAM1"]:
    dd = np.array([lb.loc[g, f"{p}_GG_counts"] - lb.loc[g, f"{p}_DN_counts"] for p in BULK7])
    bt.append(dict(cohort="GSE288479 bulk（组织层面）", gene=g, n_pairs=7, mean_diff_log2cpm=float(dd.mean()), p_one_sided_unadj=wil(dd), n_up=int((dd > 0).sum())))
BT = pd.DataFrame(bt); log("H4 在 GSE288479 bulk（未校正单侧，GG vs DN）：\n" + BT.round(4).to_string(index=False))
T.to_csv(f"{OUT}/H4_peng_stepwise.tsv", sep="\t", index=False); VT.to_csv(f"{OUT}/H4_GSE164789.tsv", sep="\t", index=False); BT.to_csv(f"{OUT}/H4_GSE288479_bulk.tsv", sep="\t", index=False)
log("===== 第 4 步完成 =====")
