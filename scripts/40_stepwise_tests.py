#!/usr/bin/env python
"""第四部分 第 1 步（分析 A）：分步检验。预设见 docs/第四部分_正常到前体驱动因素预设方案.md 第 1 节（用户已确认）。
- 两步：正常 → AAH（8 对），正常 → AIS/MIA（AIS+MIA 同患者合并）。单位 = 患者 × 组织学组；每侧 >=50 个细胞；单细胞最低 8 对，不足标「无法进行」。
- 指标（第 1.2 节，不新增）：KAC 比例、CAS 比例（AT2 谱系）、IL1B+ 巨噬细胞比例；IL1R1、CD47、RELA、RELB、NFKB1（AT2 谱系）、SIRPA（巨噬细胞）、IL1B（髓系）表达。
- 阈值直接读取第二阶段固定值（results/20_stage2/discovery/run_log.txt），不重算：KAC 0.1196、CAS 0.8261、IL1B 0.5904。
- 比例：配对 Wilcoxon 单侧（前体 > 正常）。表达：pydeseq2 `~ patient + condition`（单侧 p 由 Wald p 得到，方向=前体更高），另报告 log1p CPM 配对 Wilcoxon 单侧。
  Holm：每个指标在两步上共 2 个检验；表达指标以 pydeseq2 单侧 p 为主检验（Wilcoxon 另行 Holm 报告）。「支持」= 方向一致、单侧 p<0.05、Holm 后 p<0.05。
- AAH 不产出新候选基因；pydeseq2 全基因组结果只用于读取预设基因，不做 FDR 筛选清单（FDR 列仅为预设基因在该次运行中的 BH）。
- 「排除 II–III 期」在 AAH 组标「无法验证」；另报告去掉 P24、P25 后 6 对的方向（只报方向，不报 P）。
- ssGSEA（KAC/NF-κB，AT2 谱系伪批量；gseapy 默认参数）：只报告未校正 P，不进入 Holm 家族。
- bulk（第 1.4 节）：GSE102511 AAH vs NL（15 对）、GSE319666 AAH vs N / AIS/MIA vs N、GSE288479 GG vs DN（7 对）；
  表达类指标 + KAC/CAS ssGSEA，配对 Wilcoxon 单侧；不加 BH（报告未校正 P）。
本研究不推断干预方向的正负（表达高低不能决定抑制有益还是有害）。
输出 results/40_drivers/stepwise/；日志 results/40_drivers/run_log.txt"""
import os, re, sys, time, importlib.util, warnings
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import gseapy as gp
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/40_drivers/stepwise"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("m14", f"{P}/scripts/14_microenv_layer1.py"); m14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m14)
LOG = open(f"{P}/results/40_drivers/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()
NOTE = "本研究不推断干预方向的正负"
MIN_CELLS, MIN_PAIRS_SC = 50, 8
DROP_II_III = ["P24", "P25"]           # AAH 组中的 II–III 期患者里属 MD Anderson 者；完整排除名单 P11,P20,P24,P25 在 AAH 组仅剩 4 对（无法验证）
EXCL_AAH_STAGE = ["P11", "P20", "P24", "P25"]
LINS = {"AT2L": ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"],
        "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"], "Macrophage": ["Macrophage"]}
EXPR = [("IL1R1", "AT2L"), ("CD47", "AT2L"), ("RELA", "AT2L"), ("RELB", "AT2L"), ("NFKB1", "AT2L"), ("SIRPA", "Macrophage"), ("IL1B", "Myeloid")]
txt = open(f"{P}/results/20_stage2/discovery/run_log.txt", encoding="utf-8").read()
m = re.search(r"阈值（Normal AT2 谱系细胞第 90 百分位）：\{'KAC': np\.float64\(([0-9.]+)\), 'CAS': np\.float64\(([0-9.]+)\)\}", txt)
THR_KAC, THR_CAS = float(m.group(1)), float(m.group(2))
THR_IL1B = float(re.search(r"IL1B 第 90 百分位（Normal 巨噬细胞）= ([0-9.]+)", txt).group(1))
kac = [l.strip() for l in open(f"{P}/results/20_stage2/peng_TableS3_KAC_human.txt")]
nfkb = [l.strip() for l in open(f"{P}/results/20_stage2/peng_TableS3_NFkB_human.txt")]
CAS = ["SCGB3A2", "CEACAM6", "SPINK1", "CEACAM5"]
log(f"===== 第四部分 第 1 步 分析 A（分步检验）；固定阈值 KAC={THR_KAC:.6f} CAS={THR_CAS:.6f} IL1B={THR_IL1B} =====")

def wil(d):
    d = np.asarray(d, float); d = d[~np.isnan(d)]
    if len(d) < 5 or np.all(d == 0): return np.nan
    return float(wilcoxon(d, alternative="greater").pvalue)

a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
o = a.obs; o["patient"] = o["patient"].astype(str); o["stage"] = o["stage"].astype(str); o["fl"] = o["final_label"].astype(str)
epcam = np.asarray(a[:, "EPCAM"].layers["counts"].todense()).ravel()
X = a.layers["counts"].tocsr()
STEPS = {"AAH": ["AAH"], "AISMIA": ["AIS", "MIA"]}
patients = sorted(o["patient"].unique(), key=lambda x: int(x[1:]))

def units(labs, epcam_filter, step):
    m_ = o["fl"].isin(labs).values & ((epcam < 2) if epcam_filter else True)
    idxs, meta = {}, []
    for p in patients:
        sp_ = o["patient"].values == p
        nrm = np.where(m_ & sp_ & (o["stage"].values == "Normal"))[0]
        pre = np.where(m_ & sp_ & np.isin(o["stage"].values, STEPS[step]))[0]
        if len(nrm) < MIN_CELLS or len(pre) < MIN_CELLS: continue
        for g, ix in [("Normal", nrm), ("Precursor", pre)]:
            idxs[f"{p}_{g}"] = ix; meta.append(dict(sample=f"{p}_{g}", patient=p, condition=g, n_cells=len(ix)))
    return idxs, pd.DataFrame(meta).set_index("sample") if meta else pd.DataFrame()

# ---------- 细胞水平评分（与第二阶段相同做法；仅用固定阈值）----------
m_at = o["fl"].isin(LINS["AT2L"]).values
s2 = a[m_at].copy(); s2.obs = s2.obs[["patient", "stage"]].copy()
for nm, genes in [("KAC", kac), ("CAS", CAS)]:
    sc.tl.score_genes(s2, [g for g in genes if g in s2.var_names], score_name=nm, random_state=0)
nm_mask = (s2.obs["stage"] == "Normal").values
for nm, thr in [("KAC", THR_KAC), ("CAS", THR_CAS)]:
    rec = np.percentile(s2.obs.loc[nm_mask, nm], 90); log(f"重现性核对 {nm}：Normal P90={rec:.10f} 与固定阈值差 {abs(rec-thr):.2e}")
pos_map = -np.ones(a.n_obs, dtype=int); pos_map[np.where(m_at)[0]] = np.arange(m_at.sum())
m_mac = (o["fl"] == "Macrophage").values & (epcam < 2)
il1b_all = np.asarray(a[:, "IL1B"].X.todense()).ravel() if hasattr(a.X, "todense") else np.asarray(a[:, "IL1B"].X).ravel()
log(f"IL1B+ 比例核对：Normal 巨噬细胞 IL1B>{THR_IL1B} 比例 = {(il1b_all[m_mac & (o['stage'].values=='Normal')] > THR_IL1B).mean()*100:.1f}%（第二阶段报告 10.0%）")

frac_rows, expr_rows, dir_rows, ss_rows = [], [], [], []
tests = []
def add_test(step, indicator, kind, lineage, n_pairs, stat, p1, extra=None):
    d = dict(step=step, indicator=indicator, kind=kind, lineage=lineage, n_pairs=n_pairs, statistic=stat, p_one_sided=p1)
    if extra: d.update(extra)
    tests.append(d)

for step in STEPS:
    # ---- 比例类 ----
    idx_at, meta_at = units(LINS["AT2L"], False, step)
    idx_mac, meta_mac = units(LINS["Macrophage"], True, step)
    for ind, labs_idx, meta_, getter in [
        ("KAC_frac", idx_at, meta_at, lambda ix: float((s2.obs["KAC"].values[pos_map[ix]] > THR_KAC).mean())),
        ("CAS_frac", idx_at, meta_at, lambda ix: float((s2.obs["CAS"].values[pos_map[ix]] > THR_CAS).mean())),
        ("IL1Bpos_macrophage_frac", idx_mac, meta_mac, lambda ix: float((il1b_all[ix] > THR_IL1B).mean()))]:
        npairs = len(meta_) // 2
        if npairs < MIN_PAIRS_SC:
            add_test(step, ind, "fraction", "AT2L" if "mac" not in ind else "Macrophage", npairs, np.nan, np.nan, dict(status=f"无法进行（{npairs} 对 <{MIN_PAIRS_SC}）")); continue
        vals = pd.Series({k: getter(ix) for k, ix in labs_idx.items()})
        df = pd.DataFrame({"normal": [vals[f"{p}_Normal"] for p in meta_.patient.unique()], "precursor": [vals[f"{p}_Precursor"] for p in meta_.patient.unique()]}, index=meta_.patient.unique())
        df["diff"] = df.precursor - df.normal
        for p, r in df.iterrows(): frac_rows.append(dict(step=step, indicator=ind, patient=p, normal=r.normal, precursor=r.precursor, diff=r["diff"]))
        add_test(step, ind, "fraction", "AT2L" if "mac" not in ind else "Macrophage", npairs, float(df["diff"].mean()), wil(df["diff"]),
                 dict(n_up=int((df["diff"] > 0).sum()), status="已检验"))
        if step == "AAH":
            dd = df.drop(index=[p for p in DROP_II_III if p in df.index])["diff"]
            dir_rows.append(dict(indicator=ind, kind="fraction", n_pairs_after_drop=len(dd), mean_diff=float(dd.mean()), n_up=int((dd > 0).sum()), direction="前体更高" if dd.mean() > 0 else "前体更低或持平",
                                 note="仅方向，不报 P；去掉 P24、P25；完整排除 II–III 期（P11,P20,P24,P25）后仅 4 对：无法验证"))
    # ---- 表达类（pydeseq2 + log1p CPM Wilcoxon）----
    for lin in LINS:
        genes_here = [g for g, l in EXPR if l == lin]
        if not genes_here: continue
        idxs, meta = units(LINS[lin], lin != "AT2L", step); npairs = len(meta) // 2
        if npairs < MIN_PAIRS_SC:
            for g in genes_here: add_test(step, g, "expression", lin, npairs, np.nan, np.nan, dict(status=f"无法进行（{npairs} 对 <{MIN_PAIRS_SC}）"))
            continue
        counts = pd.DataFrame({k: np.asarray(X[ix].sum(axis=0)).ravel() for k, ix in idxs.items()}, index=a.var_names).T.astype(int)
        meta2 = meta.copy(); meta2["condition"] = pd.Categorical(meta2["condition"], categories=["Normal", "Precursor"])
        r = m14.deseq(counts, meta2[["patient", "condition"]], "~ patient + condition", ["condition", "Precursor", "Normal"])
        r["p1_up"] = np.where(r["log2FoldChange"] > 0, r["pvalue"] / 2, 1 - r["pvalue"] / 2)
        r[["log2FoldChange", "pvalue", "padj", "p1_up"]].loc[[g for g in genes_here if g in r.index]].to_csv(f"{OUT}/deseq_presetgenes_{lin}_{step}.tsv", sep="\t")
        meta.to_csv(f"{OUT}/units_{lin}_{step}.tsv", sep="\t")
        lg = np.log1p(counts.div(counts.sum(axis=1), axis=0) * 1e6)
        pts = list(meta.patient.unique())
        for g in genes_here:
            if g not in lg.columns or g not in r.index or pd.isna(r.loc[g, "pvalue"]):
                add_test(step, g, "expression", lin, npairs, np.nan, np.nan, dict(status="过滤后无检测值")); continue
            dif = pd.Series({p: lg.loc[f"{p}_Precursor", g] - lg.loc[f"{p}_Normal", g] for p in pts})
            add_test(step, g, "expression", lin, npairs, float(r.loc[g, "log2FoldChange"]), float(r.loc[g, "p1_up"]),
                     dict(p_two_sided=float(r.loc[g, "pvalue"]), padj_in_run=float(r.loc[g, "padj"]), wilcoxon_p1=wil(dif), mean_dlog1pcpm=float(dif.mean()), n_up=int((dif > 0).sum()), status="已检验"))
            if step == "AAH":
                dd = dif.drop(index=[p for p in DROP_II_III if p in dif.index])
                dir_rows.append(dict(indicator=g, kind="expression", n_pairs_after_drop=len(dd), mean_diff=float(dd.mean()), n_up=int((dd > 0).sum()), direction="前体更高" if dd.mean() > 0 else "前体更低或持平",
                                     note="仅方向（log1p CPM 差）；去掉 P24、P25"))
        # ---- ssGSEA（AT2L 伪批量）----
        if lin == "AT2L":
            ld = lg.T
            gs = {"KAC": [g for g in kac if g in ld.index], "NFKB": [g for g in nfkb if g in ld.index]}
            ss = gp.ssgsea(data=ld, gene_sets=gs, outdir=None, no_plot=True, threads=4, min_size=15)
            sres = ss.res2d.pivot(index="Name", columns="Term", values="NES").astype(float)
            for term in sres.columns:
                dd = pd.Series({p: sres.loc[f"{p}_Precursor", term] - sres.loc[f"{p}_Normal", term] for p in pts})
                ss_rows.append(dict(step=step, score=f"ssGSEA {term}", n_pairs=npairs, mean_diff=float(dd.mean()), p_one_sided_unadj=wil(dd), n_up=int((dd > 0).sum())))
    log(f"[{step}] 完成：AT2L {len(meta_at)//2} 对；巨噬细胞 {len(meta_mac)//2} 对")

T = pd.DataFrame(tests)
# Holm：每个指标在两步上共 2 个检验（主检验 p）；Wilcoxon 另行
T["p_holm_indicator"] = np.nan; T["wilcoxon_p_holm_indicator"] = np.nan
for ind, g in T.groupby("indicator"):
    ok = g["p_one_sided"].notna()
    if ok.sum(): T.loc[g.index[ok.values], "p_holm_indicator"] = multipletests(g.loc[ok, "p_one_sided"], method="holm")[1]
    if "wilcoxon_p1" in g:
        okw = g["wilcoxon_p1"].notna()
        if okw.sum(): T.loc[g.index[okw.values], "wilcoxon_p_holm_indicator"] = multipletests(g.loc[okw, "wilcoxon_p1"], method="holm")[1]
T["supported"] = (T["statistic"] > 0) & (T["p_one_sided"] < 0.05) & (T["p_holm_indicator"] < 0.05)
T["note"] = NOTE
T.to_csv(f"{OUT}/stepwise_sc_tests.tsv", sep="\t", index=False)
pd.DataFrame(frac_rows).to_csv(f"{OUT}/stepwise_fractions_per_patient.tsv", sep="\t", index=False)
pd.DataFrame(dir_rows).to_csv(f"{OUT}/AAH_drop_P24_P25_direction.tsv", sep="\t", index=False)
pd.DataFrame(ss_rows).to_csv(f"{OUT}/stepwise_ssgsea_sc.tsv", sep="\t", index=False)
pd.set_option("display.width", 260, "display.max_colwidth", 40)
log("单细胞分步检验（" + NOTE + "）：\n" + T[["step", "indicator", "kind", "n_pairs", "statistic", "p_one_sided", "p_holm_indicator", "supported", "wilcoxon_p1", "n_up", "status"]].round(5).to_string(index=False))
log("AAH 去掉 P24、P25（只报方向）：\n" + pd.DataFrame(dir_rows)[["indicator", "n_pairs_after_drop", "mean_diff", "n_up", "direction"]].round(4).to_string(index=False))
log("ssGSEA（单细胞伪批量，未校正，不进 Holm）：\n" + pd.DataFrame(ss_rows).round(4).to_string(index=False))
del a, X, s2

# ====================== bulk ======================
GENES_B = [g for g, _ in EXPR]
def bulk_tests(A, B, name, cohort):
    """A、B：基因×患者，同列顺序；单侧 A>B 配对 Wilcoxon"""
    rows = []
    for g in GENES_B:
        if g not in A.index: rows.append(dict(cohort=cohort, comparison=name, indicator=g, kind="expression", n_pairs=A.shape[1], status="无此基因")); continue
        d = A.loc[g].values.astype(float) - B.loc[g].values.astype(float)
        rows.append(dict(cohort=cohort, comparison=name, indicator=g, kind="expression", n_pairs=len(d), mean_diff=float(d.mean()), p_one_sided_unadj=wil(d), n_up=int((d > 0).sum()), status="已检验"))
    gs = {"KAC": [g for g in kac if g in A.index], "CAS": [g for g in CAS if g in A.index]}
    both = pd.concat([A.add_suffix("|P"), B.add_suffix("|N")], axis=1)
    ss = gp.ssgsea(data=both, gene_sets=gs, outdir=None, no_plot=True, threads=4, min_size=3)
    sres = ss.res2d.pivot(index="Name", columns="Term", values="NES").astype(float)
    for term in sres.columns:
        d = np.array([sres.loc[f"{c}|P", term] - sres.loc[f"{c}|N", term] for c in A.columns])
        rows.append(dict(cohort=cohort, comparison=name, indicator=f"ssGSEA {term}（{len(gs[term])} 个基因）", kind="ssGSEA", n_pairs=len(d), mean_diff=float(d.mean()), p_one_sided_unadj=wil(d), n_up=int((d > 0).sum()), status="已检验"))
    return rows
brows = []
# GSE102511（文件已是 log2 尺度表达）
b1 = pd.read_csv(f"{P}/raw/bulk/GSE102511_Smruthy_etal_allsamples_TPM.txt.gz", sep="\t", index_col=0, compression="gzip").groupby(level=0).mean()
don = sorted({c.split("-")[0] for c in b1.columns}); pairs1 = [d_ for d_ in don if f"{d_}-NL" in b1.columns and f"{d_}-AAH" in b1.columns]
A1 = b1[[f"{d_}-AAH" for d_ in pairs1]]; B1 = b1[[f"{d_}-NL" for d_ in pairs1]]; A1.columns = B1.columns = pairs1
brows += bulk_tests(A1, B1, "AAH vs NL", "GSE102511"); log(f"GSE102511：{len(pairs1)} 对")
# GSE319666
sm = pd.read_csv(f"{P}/results/20_stage2/gse319666_sample_meta.tsv", sep="\t")
cnt = pd.read_csv(f"{P}/raw/bulk_GSE319666_HTAN_PML/GSE319666_counts.tsv.gz", sep="\t", index_col=0, compression="gzip")
gtf = pd.read_csv(f"{P}/raw/annotation/gencode.v44.basic.annotation.gtf.gz", sep="\t", comment="#", header=None, usecols=[2, 8], names=["f", "a"], low_memory=False)
gtf = gtf[gtf["f"] == "gene"]; mp = pd.Series(gtf["a"].str.extract(r'gene_name "([^"]+)"')[0].values, index=gtf["a"].str.extract(r'gene_id "([^".]+)')[0].values); mp = mp[~mp.index.duplicated()]
cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True); cnt["g"] = mp.reindex(cnt.index).values; cnt = cnt[cnt["g"].notna()].groupby("g").sum()
lc = np.log2(cnt.div(cnt.sum(axis=0), axis=1) * 1e6 + 1); sm = sm[sm["col"].isin(lc.columns)]
def pmean(h):
    s = sm[sm["hist"] == h]; return lc[s["col"].values].T.groupby(s["patient"].values).mean().T
Ng, AAH, AM = pmean("N"), pmean("AAH"), pmean("AIS/MIA")
pa = [p for p in AAH.columns if p in Ng.columns]; pm = [p for p in AM.columns if p in Ng.columns]
brows += bulk_tests(AAH[pa], Ng[pa], "AAH vs N", "GSE319666"); brows += bulk_tests(AM[pm], Ng[pm], "AIS/MIA vs N", "GSE319666")
log(f"GSE319666：AAH vs N {len(pa)} 对；AIS/MIA vs N {len(pm)} 对")
# GSE288479 bulk（7 对 GG vs DN；组织层面方向参考）
BULK7 = ["SMC06", "SMC08", "SMC11", "SMC13", "SMC18", "SMC19", "SMC20"]
bk = pd.read_csv(f"{P}/raw/multi_GSE288479_partsolid/GSE288479_raw_counts_All_samples.txt.gz", index_col=0, compression="gzip")
bk = bk.drop(columns=["gene_id", "gene_type"]).groupby("gene_name").sum()
bk = bk[[f"{p}_{s}_counts" for p in BULK7 for s in ["DN", "GG"]]]
lb = np.log2(bk.div(bk.sum(axis=0), axis=1) * 1e6 + 1)
Ab = lb[[f"{p}_GG_counts" for p in BULK7]]; Bb = lb[[f"{p}_DN_counts" for p in BULK7]]; Ab.columns = Bb.columns = BULK7
brows += bulk_tests(Ab, Bb, "GG vs DN（组织层面方向参考）", "GSE288479"); log("GSE288479 bulk：7 对")
BT = pd.DataFrame(brows); BT["note"] = NOTE + "；bulk 为混合细胞，方向一致可能来自细胞构成变化；不加 BH"
BT.to_csv(f"{OUT}/stepwise_bulk_tests.tsv", sep="\t", index=False)
log("bulk 分步补充（未校正单侧 P）：\n" + BT[["cohort", "comparison", "indicator", "n_pairs", "mean_diff", "p_one_sided_unadj", "n_up"]].round(4).to_string(index=False))
log("===== 第 1 步完成 =====")
