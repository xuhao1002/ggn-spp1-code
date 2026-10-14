#!/usr/bin/env python
"""
第二阶段：外部验证 + bulk 临床层 + 第 5 层 + 汇总（预设见 docs/第二阶段预设方案.md，2026-10-08 用户确认，运行前写定）。
- 外部验证（任一通过即通过）：
  GSE164789：AIS/MIA 病灶（-T）vs 同患者癌旁（-A），单位=患者×组织学，pydeseq2 `~ pair + tissue`；单位两侧>=50 个细胞，可配对数<8 的细胞群「无法验证」。
  GSE288479 bulk：GG vs 同患者 DN（只用 SMC06、08、11、13、18、19、20 七名患者），pydeseq2 `~ patient + site`，>=7 对；bulk 只代表组织层面方向。
  规则：对候选基因单侧检验（病灶/GG > 对照），BH 校正（GSE164789 在细胞群内；bulk 在全部候选基因中），q<0.05 且方向一致 = 通过。
- 临床相关（替代第 4 层）：GSE102511（AAH vs NL，15 对）、GSE319666（AAH vs N，AIS/MIA vs N；同患者多个样本取平均 log-CPM）；配对 Wilcoxon 单侧（前体 > 正常），
  任一比较方向一致且 P<0.05 即支持；不加 BH（同时报告 q）。
- 第 5 层：沿用第一阶段修正版（cleavage 去掉；A 级须 ChEMBL/TTD 抑制类记录；免疫来源靶点的血液/骨髓规则）。
- 成功标准：同时满足 发现 FDR<0.05 且前体>正常、通过来源过滤与非增殖、外部验证通过、bulk 支持、第 5 层通过。
输出 results/20_stage2/validation、layers
"""
import os, sys, re, time, zipfile, importlib.util, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import wilcoxon, spearmanr
from statsmodels.stats.multitest import multipletests
import gseapy as gp
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); R = f"{P}/results/20_stage2"
OUTV, OUTL = f"{R}/validation", f"{R}/layers"; os.makedirs(OUTV, exist_ok=True); os.makedirs(OUTL, exist_ok=True)
sys.path.insert(0, f"{P}/scripts")
spec = importlib.util.spec_from_file_location("m14", f"{P}/scripts/14_microenv_layer1.py"); m14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m14)
LOG = open(f"{R}/stage2_run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()
VL = {"AT2L": ["AT2", "Club/secretory"], "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"], "Macrophage": ["Macrophage"],
      "Fibroblast": ["Fibroblast"], "T/NK": ["T cell", "NK"], "Endothelial": ["Endothelial"]}
EPI_V = ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal"]
MIN_PAIRS_SC, MIN_PAIRS_BULK = 8, 7
BULK7 = ["SMC06", "SMC08", "SMC11", "SMC13", "SMC18", "SMC19", "SMC20"]
kac = [l.strip() for l in open(f"{R}/peng_TableS3_KAC_human.txt")]; nfkb = [l.strip() for l in open(f"{R}/peng_TableS3_NFkB_human.txt")]
CAS = ["SCGB3A2", "CEACAM6", "SPINK1", "CEACAM5"]
cand = pd.read_csv(f"{R}/discovery/candidates_stage2_eligible.tsv", sep="\t", index_col=0)
cand["gene"] = cand.index
log(f"合格候选（基因×细胞群）{len(cand)}：{cand['lineage'].value_counts().to_dict()}")

def one_sided(r, positive=True):
    return np.where(r["log2FoldChange"] > 0, r["pvalue"] / 2, 1 - r["pvalue"] / 2)

# ====================== GSE164789 ======================
d = f"{P}/results/05_validation/GSE164789"
v = sc.read_h5ad(f"{d}/GSE164789_annotated.h5ad")
ann = pd.read_csv(f"{d}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
v.obs["fl"] = v.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
for c in ["patient", "tissue", "histology", "sample"]: v.obs[c] = v.obs[c].astype(str)
ep = np.asarray(v[:, "EPCAM"].layers["counts"].todense()).ravel(); pt = np.asarray(v[:, "PTPRC"].layers["counts"].todense()).ravel()
isepi = v.obs["fl"].isin(EPI_V).values
clean = ~(isepi & (pt >= 2)) & ~(~isepi & (ep >= 2))
v = v[clean & v.obs["histology"].isin(["AIS", "MIA"]).values].copy()
v.obs["unit"] = v.obs["patient"] + "_" + v.obs["histology"]
log(f"GSE164789：AIS/MIA 细胞 {v.n_obs}；单位 {v.obs['unit'].nunique()}")
Xv = v.layers["counts"].tocsr()
def v_units(labs):
    m = v.obs["fl"].isin(labs).values; idxs, meta = {}, []
    for u in sorted(v.obs["unit"].unique()):
        a_ = np.where(m & (v.obs["unit"].values == u) & (v.obs["tissue"].values == "A"))[0]; t_ = np.where(m & (v.obs["unit"].values == u) & (v.obs["tissue"].values == "T"))[0]
        if len(a_) < 50 or len(t_) < 50: continue
        for g, ix in [("A", a_), ("T", t_)]: idxs[f"{u}|{g}"] = ix; meta.append(dict(sample=f"{u}|{g}", pair=u, tissue=g, n_cells=len(ix)))
    return idxs, pd.DataFrame(meta).set_index("sample") if meta else pd.DataFrame()
val_res, val_status = {}, {}
for lin in VL:
    idxs, meta = v_units(VL[lin]); npairs = len(meta) // 2
    val_status[lin] = npairs
    if npairs < 2: log(f"[GSE164789][{lin}] 可配对 {npairs}"); continue
    counts = pd.DataFrame({k: np.asarray(Xv[ix].sum(axis=0)).ravel() for k, ix in idxs.items()}, index=v.var_names).T.astype(int)
    meta["tissue"] = pd.Categorical(meta["tissue"], categories=["A", "T"])
    r = m14.deseq(counts, meta[["pair", "tissue"]], "~ pair + tissue", ["tissue", "T", "A"]); r["p1_up"] = one_sided(r)
    r.to_csv(f"{OUTV}/GSE164789_deg_{lin.replace('/','_')}.tsv", sep="\t"); meta.to_csv(f"{OUTV}/GSE164789_units_{lin.replace('/','_')}.tsv", sep="\t")
    val_res[lin] = r; log(f"[GSE164789][{lin}] 可配对 {npairs} 对；检测基因 {int(r['padj'].notna().sum())}" + ("" if npairs >= MIN_PAIRS_SC else "；<8 对：无法验证"))

# --- H1/H2/H3 在 GSE164789 ---
def hyp_gse164789():
    res = {}
    tests = []
    def deseq_row(lin, gene):
        r = val_res.get(lin)
        if r is None or gene not in r.index or pd.isna(r.loc[gene, "pvalue"]): return (np.nan,) * 3
        return r.loc[gene, "log2FoldChange"], r.loc[gene, "padj"], r.loc[gene, "p1_up"]
    for fam, tid, lin, gene, prim in [("H1", "H1a", "AT2L", "IL1R1", True), ("H1", "H1b", "Myeloid", "IL1B", False), ("H1", "H1c-RELA", "AT2L", "RELA", False), ("H1", "H1c-RELB", "AT2L", "RELB", False),
                                      ("H1", "H1c-NFKB1", "AT2L", "NFKB1", False), ("H2", "H2a", "AT2L", "CD47", False), ("H2", "H2b", "Macrophage", "SIRPA", False)]:
        l, q, p1 = deseq_row(lin, gene); tests.append(dict(family=fam, test=tid, primary=prim, description=f"{lin} {gene} 伪批量", statistic=l, n=val_status.get(lin, 0), p_one_sided=p1, note=f"双侧 FDR={q:.3g}" if pd.notna(q) else ""))
    # 细胞比例（阈值取对照组织 -A 的细胞）
    m_at = v.obs["fl"].isin(VL["AT2L"]).values; s2 = v[m_at].copy()
    for nm, genes in [("KAC", kac), ("CAS", CAS)]:
        sc.tl.score_genes(s2, [g for g in genes if g in s2.var_names], score_name=nm, random_state=0)
    a_mask = (s2.obs["tissue"] == "A").values
    thr = {nm: np.percentile(s2.obs.loc[a_mask, nm], 90) for nm in ["KAC", "CAS"]}
    mm = v[(v.obs["fl"] == "Macrophage").values]
    il1b = np.asarray(mm[:, "IL1B"].X.todense()).ravel() if sp.issparse(mm.X) else np.asarray(mm[:, "IL1B"].X).ravel(); sirpa = np.asarray(mm[:, "SIRPA"].layers["counts"].todense()).ravel()
    thr_i = np.percentile(il1b[(mm.obs["tissue"] == "A").values], 90); pos = (il1b > thr_i) if thr_i > 0 else (il1b > 0)
    rows = []
    for u in sorted(v.obs["unit"].unique()):
        for g in ["A", "T"]:
            sa = np.where((s2.obs["unit"].values == u) & (s2.obs["tissue"].values == g))[0]; sm = np.where((mm.obs["unit"].values == u) & (mm.obs["tissue"].values == g))[0]
            rows.append(dict(unit=u, tissue=g, n_at2=len(sa), n_mac=len(sm), KAC=float((s2.obs["KAC"].values[sa] > thr["KAC"]).mean()) if len(sa) >= 50 else np.nan,
                             CAS=float((s2.obs["CAS"].values[sa] > thr["CAS"]).mean()) if len(sa) >= 50 else np.nan,
                             IL1B=float(pos[sm].mean()) if len(sm) >= 50 else np.nan, SIRPA=float((sirpa[sm] >= 1).mean()) if len(sm) >= 50 else np.nan))
    fr = pd.DataFrame(rows); fr.to_csv(f"{OUTV}/GSE164789_cell_fractions.tsv", sep="\t", index=False)
    def paired(col, desc, tid, fam, prim):
        t = fr[fr.tissue == "T"].set_index("unit")[col]; a = fr[fr.tissue == "A"].set_index("unit")[col]; dd = (t - a.reindex(t.index)).dropna()
        tests.append(dict(family=fam, test=tid, primary=prim, description=desc, statistic=float(dd.mean()) if len(dd) else np.nan, n=len(dd),
                          p_one_sided=wilcoxon(dd, alternative="greater").pvalue if len(dd) >= 5 else np.nan, note=f"病灶更高 {int((dd>0).sum())}/{len(dd)}；阈值取 -A 的第 90 百分位（本研究自定义）"))
    paired("KAC", "KAC 细胞比例（AT2 谱系）", "H1d", "H1", False); paired("IL1B", "IL1B+ 巨噬细胞比例", "H1e", "H1", True)
    paired("SIRPA", "SIRPA+ 巨噬细胞比例", "H2c", "H2", False); paired("CAS", "CAS 样上皮比例（AT2 谱系）", "H3", "H3", False)
    T = pd.DataFrame(tests); T["p_holm_family"] = np.nan
    for fam, g in T.groupby("family"):
        ok = g["p_one_sided"].notna()
        if ok.sum(): T.loc[g.index[ok.values], "p_holm_family"] = multipletests(g.loc[ok, "p_one_sided"], method="holm")[1]
    T["supported"] = (T["statistic"] > 0) & (T["p_one_sided"] < 0.05) & (T["p_holm_family"] < 0.05)
    T.to_csv(f"{OUTV}/GSE164789_hypotheses.tsv", sep="\t", index=False)
    log("GSE164789 假说检验：\n" + T[["family", "test", "description", "statistic", "n", "p_one_sided", "p_holm_family", "supported"]].round(5).to_string(index=False))
hyp_gse164789()
del v, Xv

# --- 候选验证 ---
def sc_validate(c):
    out = []
    for lin, g in c.groupby("lineage"):
        r = val_res.get(lin); npairs = val_status.get(lin, 0)
        g = g.copy()
        if r is None or npairs < MIN_PAIRS_SC:
            g["val164789_status"] = f"无法验证（可配对 {npairs} 对 <{MIN_PAIRS_SC}）"; g["val164789_pass"] = False
        else:
            rr = r.reindex(g.index); ok = rr["pvalue"].notna().values
            p1 = rr["p1_up"].values; q = np.full(len(g), np.nan)
            if ok.sum(): q[ok] = multipletests(p1[ok], method="fdr_bh")[1]
            g["val164789_lfc"] = rr["log2FoldChange"].values; g["val164789_p1"] = p1; g["val164789_q"] = q
            g["val164789_pass"] = (rr["log2FoldChange"].values > 0) & (q < 0.05)
            g["val164789_status"] = np.where(ok, np.where(g["val164789_pass"], "通过", "未通过"), "无检测值")
        out.append(g)
    return pd.concat(out)
cand = sc_validate(cand)
for lin, g in cand.groupby("lineage"):
    log(f"[GSE164789 验证][{lin}] 候选 {len(g)}；{g['val164789_status'].value_counts().to_dict()}")

# ====================== GSE288479 bulk ======================
bk = pd.read_csv(f"{P}/raw/multi_GSE288479_partsolid/GSE288479_raw_counts_All_samples.txt.gz", index_col=0, compression="gzip")
bk = bk.drop(columns=["gene_id", "gene_type"]).groupby("gene_name").sum()
cols = {f"{p}_{s}_counts": (p, s) for p in BULK7 for s in ["DN", "GG", "S"]}
bk = bk[list(cols)]; bk.columns = [f"{p}_{s}" for p, s in cols.values()]
def bulk_de(sites, ref, name):
    sel = [c for c in bk.columns if c.split("_")[1] in sites]
    meta = pd.DataFrame({"patient": [c.split("_")[0] for c in sel], "site": [c.split("_")[1] for c in sel]}, index=sel)
    meta["site"] = pd.Categorical(meta["site"], categories=[ref, [s for s in sites if s != ref][0]])
    r = m14.deseq(bk[sel].T.astype(int), meta, "~ patient + site", ["site", [s for s in sites if s != ref][0], ref]); r["p1_up"] = one_sided(r)
    r.to_csv(f"{OUTV}/GSE288479_bulk_{name}.tsv", sep="\t"); return r
gg = bulk_de(["DN", "GG"], "DN", "GG_vs_DN"); sg = bulk_de(["GG", "S"], "GG", "S_vs_GG")
log(f"GSE288479 bulk：7 对 GG vs DN；检测基因 {int(gg['padj'].notna().sum())}；FDR<0.05 {int((gg['padj']<0.05).sum())}（GG 更高 {int(((gg['padj']<0.05)&(gg.log2FoldChange>0)).sum())}）；S vs GG FDR<0.05 {int((sg['padj']<0.05).sum())}")
genes_u = sorted(set(cand.index)); rr = gg.reindex(genes_u); ok = rr["pvalue"].notna().values
q288 = pd.Series(np.nan, index=genes_u)
if ok.sum(): q288[ok] = multipletests(rr["p1_up"].values[ok], method="fdr_bh")[1]
cand["val288479_lfc"] = gg["log2FoldChange"].reindex(cand.index).values; cand["val288479_p1"] = gg["p1_up"].reindex(cand.index).values
cand["val288479_q"] = q288.reindex(cand.index).values
cand["val288479_pass"] = (cand["val288479_lfc"] > 0) & (cand["val288479_q"] < 0.05)
cand["val288479_status"] = np.where(cand["val288479_lfc"].notna(), np.where(cand["val288479_pass"], "通过", "未通过"), "无检测值")
cand["valS_vs_GG_lfc"] = sg["log2FoldChange"].reindex(cand.index).values
log(f"[GSE288479 bulk 验证] 参与校正基因 {int(ok.sum())}（候选唯一基因 {len(genes_u)}）；通过 {int(cand['val288479_pass'].sum())}（基因×群行）")
cand["external_validation_pass"] = cand["val164789_pass"] | cand["val288479_pass"]
# bulk 假说表达（只报告）
rep = gg.reindex(["IL1R1", "IL1B", "RELA", "RELB", "NFKB1", "CD47", "SIRPA", "KRT8"])[["log2FoldChange", "pvalue", "padj", "p1_up"]]
rep.to_csv(f"{OUTV}/GSE288479_bulk_hypothesis_genes.tsv", sep="\t"); log("GSE288479 bulk GG vs DN（组织水平）\n" + rep.round(4).to_string())

# ====================== bulk 临床层 ======================
b1 = pd.read_csv(f"{P}/raw/bulk/GSE102511_Smruthy_etal_allsamples_TPM.txt.gz", sep="\t", index_col=0, compression="gzip")
b1 = b1.groupby(level=0).mean()
don = sorted({c.split("-")[0] for c in b1.columns}); pairs1 = [d_ for d_ in don if f"{d_}-NL" in b1.columns and f"{d_}-AAH" in b1.columns]
def paired_test(A, B, genes):
    """A、B：基因×患者（同列顺序），单侧 A>B 的配对 Wilcoxon"""
    out = {}
    for g in genes:
        if g not in A.index: continue
        d_ = A.loc[g].values.astype(float) - B.loc[g].values.astype(float)     # 2026-10-08 修正：按位置相减（两侧列名不同，按标签对齐会得到 NaN）
        if np.all(d_ == 0): out[g] = (0.0, 1.0); continue
        out[g] = (float(d_.mean()), wilcoxon(d_, alternative="greater").pvalue)
    return pd.DataFrame(out, index=["mean_diff", "p1"]).T
genes_all = sorted(set(cand.index))
A1 = b1[[f"{d_}-AAH" for d_ in pairs1]]; B1 = b1[[f"{d_}-NL" for d_ in pairs1]]
t1 = paired_test(A1, B1, genes_all); log(f"GSE102511：{len(pairs1)} 对（AAH vs NL）；值为文件给出的 log2 尺度表达")
sm = pd.read_csv(f"{R}/gse319666_sample_meta.tsv", sep="\t")
cnt = pd.read_csv(f"{P}/raw/bulk_GSE319666_HTAN_PML/GSE319666_counts.tsv.gz", sep="\t", index_col=0, compression="gzip")
gtf = pd.read_csv(f"{P}/raw/annotation/gencode.v44.basic.annotation.gtf.gz", sep="\t", comment="#", header=None, usecols=[2, 8], names=["f", "a"], low_memory=False)
gtf = gtf[gtf["f"] == "gene"]; mp = pd.Series(gtf["a"].str.extract(r'gene_name "([^"]+)"')[0].values, index=gtf["a"].str.extract(r'gene_id "([^".]+)')[0].values); mp = mp[~mp.index.duplicated()]
cnt.index = cnt.index.str.replace(r"\.\d+$", "", regex=True); cnt["g"] = mp.reindex(cnt.index).values; cnt = cnt[cnt["g"].notna()].groupby("g").sum()
lc = np.log2(cnt.div(cnt.sum(axis=0), axis=1) * 1e6 + 1)
sm = sm[sm["col"].isin(lc.columns)]      # 2026-10-08 修正：计数文件列名取自 GEO 样本标题方括号内的 ID（119 个样本的 htan id 字段与列名前缀不同）
def patient_mean(h):
    s = sm[sm["hist"] == h]; return lc[s["col"].values].T.groupby(s["patient"].values).mean().T
Ngr, AAH, AM = patient_mean("N"), patient_mean("AAH"), patient_mean("AIS/MIA")
pa = [p for p in AAH.columns if p in Ngr.columns]; pm = [p for p in AM.columns if p in Ngr.columns]
t2 = paired_test(AAH[pa], Ngr[pa], genes_all); t3 = paired_test(AM[pm], Ngr[pm], genes_all)
log(f"GSE319666：AAH vs N {len(pa)} 对；AIS/MIA vs N {len(pm)} 对（同患者多个样本先取平均 log2 CPM）")
bulk = pd.DataFrame(index=genes_all)
for nm, t in [("GSE102511_AAHvsNL", t1), ("GSE319666_AAHvsN", t2), ("GSE319666_AISMIAvsN", t3)]:
    bulk[f"{nm}_diff"] = t["mean_diff"]; bulk[f"{nm}_p1"] = t["p1"]
    ok = t["p1"].notna(); qq = pd.Series(np.nan, index=t.index); qq[ok] = multipletests(t.loc[ok, "p1"], method="fdr_bh")[1]; bulk[f"{nm}_q"] = qq
    bulk[f"{nm}_support"] = (bulk[f"{nm}_diff"] > 0) & (bulk[f"{nm}_p1"] < 0.05)
bulk["bulk_support"] = bulk[[c for c in bulk.columns if c.endswith("_support")]].any(axis=1)
bulk["bulk_which"] = bulk.apply(lambda r: ",".join(c.replace("_support", "") for c in bulk.columns if c.endswith("_support") and r[c]), axis=1)
bulk.to_csv(f"{OUTV}/bulk_layer4_candidates.tsv", sep="\t")
cand = cand.join(bulk[["bulk_support", "bulk_which"]], how="left"); cand["bulk_support"] = cand["bulk_support"].fillna(False)
log(f"bulk 临床层：候选唯一基因 {len(genes_all)}；支持 {int(bulk.bulk_support.sum())}（GSE102511 {int(bulk['GSE102511_AAHvsNL_support'].sum())}，GSE319666 AAH {int(bulk['GSE319666_AAHvsN_support'].sum())}，AIS/MIA {int(bulk['GSE319666_AISMIAvsN_support'].sum())}）")

# ====================== 第 5 层 ======================
S = f"{P}/raw/safety_druggability"
dg = pd.read_csv(f"{S}/dgidb_interactions.tsv", sep="\t", low_memory=False); dg["itype"] = dg["interaction_type"].fillna("").str.lower()
SRC = ["ChEMBL", "TTD", "GuideToPharmacology", "DTC"]
INHIB = {"inhibitor", "antagonist", "blocker", "antibody", "negative modulator", "inverse agonist", "suppressor", "antisense oligonucleotide"}
UNSPEC = {"", "modulator", "binder", "other/unknown"}
dg = dg[dg["drug_name"].notna() | dg["drug_claim_name"].notna()].copy(); dg["drug"] = dg["drug_name"].fillna(dg["drug_claim_name"]); dg["appr"] = dg["approved"].astype(str).str.lower().eq("true")
inh = dg[dg["itype"].isin(INHIB)]; uns = dg[dg["itype"].isin(UNSPEC) & dg["interaction_source_db_name"].isin(SRC)]
gt = pd.read_csv(f"{S}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
core = ["Heart - Left Ventricle", "Heart - Atrial Appendage", "Liver", "Kidney - Cortex", "Kidney - Medulla"]
gt = gt.groupby("Description")[core + ["Whole Blood"]].max()
z = zipfile.ZipFile(f"{S}/normal_ihc_data.tsv.zip"); ih = pd.read_csv(z.open(z.namelist()[0]), sep="\t"); ih = ih[ih["Reliability"].isin(["Approved", "Supported", "Enhanced"])]
ctc = {"Heart muscle": ["cardiomyocytes"], "Liver": ["hepatocytes"], "Kidney": ["cells in glomeruli", "cells in tubules", "proximal tubules (cell body)", "proximal tubules (microvilli)", "distal tubules", "collecting ducts"]}
def hpa_high(tc): x = ih[ih.apply(lambda r: r["Cell type"] in tc.get(r["Tissue"], []), axis=1)]; return set(x.loc[x["Level"] == "High", "Gene name"])
hp_core = hpa_high(ctc); hp_bm = hpa_high({"Bone marrow": ["hematopoietic cells"]})
ratio_p = pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv", sep="\t", index_col=0)
imm_cols = {"Myeloid": "med_Myeloid", "T/NK": "med_T/NK"}; other_imm = ["med_T/NK", "med_Myeloid", "med_B/Plasma", "med_Mast"]
l5 = []
for (gene, lin), r in cand.groupby([cand.index, "lineage"]):
    a_ = inh[inh["gene_name"] == gene]; u_ = uns[uns["gene_name"] == gene]; allr = pd.concat([a_, u_])
    a_ok = a_[a_["appr"] & a_["interaction_source_db_name"].isin(["ChEMBL", "TTD"])]; appr = a_ok["drug"].str.upper().unique().tolist()
    oth = allr[~allr["drug"].str.upper().isin(appr)]["drug"].str.upper().unique().tolist()
    row = dict(gene=gene, lineage=lin, drug_grade="无" if allr.empty else ("A" if appr else "B"), drugs=";".join((appr + oth)[:6]), type_to_verify=(not allr.empty) and a_.empty)
    tp = gt.loc[gene] if gene in gt.index else None
    row["gtex_core_max"] = float(tp[core].max()) if tp is not None else np.nan; row["gtex_blood"] = float(tp["Whole Blood"]) if tp is not None else np.nan
    veto_core = (row["gtex_core_max"] >= 50) or (gene in hp_core)
    if lin in imm_cols:
        others = [c for c in other_imm if c != imm_cols[lin]]
        mL = ratio_p.loc[gene, imm_cols[lin]] if gene in ratio_p.index else np.nan; mo = ratio_p.loc[gene, others].max() if gene in ratio_p.index else np.nan
        fold = (mL / mo) if (mo and mo > 0) else (np.inf if mL and mL > 0 else np.nan); row["fold_vs_other_immune"] = fold
        bm_ok = (fold >= 3) and (gene not in hp_bm); veto = veto_core or not bm_ok
    else:
        bm_veto = (row["gtex_blood"] >= 50) or (gene in hp_bm); veto = veto_core or bm_veto
    row["safety_veto"] = bool(veto); row["L5_pass"] = (row["drug_grade"] in ("A", "B")) and not veto
    l5.append(row)
l5 = pd.DataFrame(l5).set_index(["gene", "lineage"])
cand = cand.reset_index(drop=True).merge(l5.reset_index(), on=["gene", "lineage"], how="left")

# ====================== 汇总与成功标准 ======================
cand["discovery_ok"] = True      # 合格候选已满足：前体>正常、FDR<0.05、通过来源过滤、非增殖
cand["stage2_pass_all"] = cand["discovery_ok"] & cand["external_validation_pass"] & cand["bulk_support"] & cand["L5_pass"]
cand.to_csv(f"{R}/stage2_candidates_integrated.tsv", sep="\t", index=False)
log("===== 汇总 =====")
for lin, g in cand.groupby("lineage"):
    log(f"[{lin}] 合格候选 {len(g)}；外部验证通过 {int(g.external_validation_pass.sum())}（GSE164789 {int(g.val164789_pass.sum())}、GSE288479 bulk {int(g.val288479_pass.sum())}）；bulk 支持 {int(g.bulk_support.sum())}；第 5 层通过 {int(g.L5_pass.sum())}；"
        f"验证&bulk {int((g.external_validation_pass & g.bulk_support).sum())}；验证&bulk&第5层 {int(g.stage2_pass_all.sum())}")
fin = cand[cand.stage2_pass_all]
log(f"第二阶段成功标准（五项同时满足）：{len(fin)} 个基因×群；基因：{sorted(set(fin.gene))}")
for _, r in fin.iterrows():
    log(f"  {r.gene}（{r.lineage}）log2FC={r.log2FoldChange:+.2f} FDR={r.padj:.2g}；GSE164789 {r.val164789_status}；GSE288479 {r.val288479_status}；bulk {r.bulk_which}；药物 {r.drug_grade} {r.drugs}")
# 通过前 4 项而第 5 层失败的（供参考）
near = cand[cand.external_validation_pass & cand.bulk_support & ~cand.L5_pass]
log(f"通过发现+验证+bulk、但第 5 层未通过：{len(near)} 个基因×群；前 25：" + ", ".join(f"{r.gene}({r.lineage},{r.drug_grade},veto={r.safety_veto})" for _, r in near.sort_values('padj').head(25).iterrows()))
log("===== 完成 =====")
