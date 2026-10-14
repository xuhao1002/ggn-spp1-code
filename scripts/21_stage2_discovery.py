#!/usr/bin/env python
"""第二阶段发现分析 + 假说检验（预设见 docs/第二阶段预设方案.md，2026-10-08 用户确认，运行前写定）。
- Peng 单核数据：前体（同患者 AAH/AIS/MIA 合并）vs 同患者 Normal（Normal+Normal1 合并）配对伪批量，pydeseq2 `~ patient + condition`，FDR<0.05；每侧>=50 个细胞。
- 细胞群：AT2 谱系（AT2 + AT2/Club-like + 增殖上皮）、髓系、成纤维、T/NK、内皮（血管）；另有巨噬细胞（Macrophage）只用于 H1e、H2。非上皮群剔除 EPCAM 计数>=2 的细胞。
- 候选：前体>正常（log2FC>0）且 FDR<0.05；不属于第二阶段增殖基因集（196 个）；通过来源过滤（上皮：第一阶段规则；微环境：反向规则）。
- 敏感性：前体仅 AIS/MIA；排除病理 II–III 期患者（P5、P11、P20、P24、P25）。只报告。
- 假说 H1、H2、H3：只报告，不入选。按假说家族 Holm 校正，同时报告未校正 P。
输出 results/20_stage2/discovery/"""
import os, sys, time, importlib.util, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import wilcoxon, spearmanr
from statsmodels.stats.multitest import multipletests
import gseapy as gp
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/20_stage2/discovery"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{P}/scripts")
import prolif_genes
spec = importlib.util.spec_from_file_location("m14", f"{P}/scripts/14_microenv_layer1.py"); m14 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m14)
LINS = {"AT2L": ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"],
        "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"],
        "Fibroblast": ["Fibroblast"], "T/NK": ["T cell", "NK", "Proliferating T/NK"], "Endothelial": ["Endothelial"],
        "Macrophage": ["Macrophage"]}
EXCL = ["P5", "P11", "P20", "P24", "P25"]
MIN_CELLS = 50
LOG = open(f"{OUT}/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()
prolif2 = prolif_genes.prolif_set_stage2()
kac = [l.strip() for l in open(f"{P}/results/20_stage2/peng_TableS3_KAC_human.txt")]
nfkb = [l.strip() for l in open(f"{P}/results/20_stage2/peng_TableS3_NFkB_human.txt")]
CAS = ["SCGB3A2", "CEACAM6", "SPINK1", "CEACAM5"]
def p1_of(r):
    return np.where(r["log2FoldChange"] > 0, r["pvalue"] / 2, 1 - r["pvalue"] / 2)

def main():
    log(f"===== 开始；第二阶段增殖基因集 {len(prolif2)} 个 =====")
    a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
    o = a.obs; o["patient"] = o["patient"].astype(str); o["stage"] = o["stage"].astype(str); o["fl"] = o["final_label"].astype(str)
    epcam = np.asarray(a[:, "EPCAM"].layers["counts"].todense()).ravel()
    # ---------- 伪批量 ----------
    def units(labs, epcam_filter, variant):
        m = o["fl"].isin(labs).values & ((epcam < 2) if epcam_filter else True)
        idxs, meta = {}, []
        for p in sorted(o["patient"].unique(), key=lambda x: int(x[1:])):
            if variant == "stageI" and p in EXCL: continue
            sp_ = o["patient"].values == p
            nrm = np.where(m & sp_ & (o["stage"].values == "Normal"))[0]
            st = ["AIS", "MIA"] if variant == "aismia" else ["AAH", "AIS", "MIA"]
            pre = np.where(m & sp_ & np.isin(o["stage"].values, st))[0]
            if len(nrm) < MIN_CELLS or len(pre) < MIN_CELLS: continue
            for g, ix in [("Normal", nrm), ("Precursor", pre)]:
                idxs[f"{p}_{g}"] = ix; meta.append(dict(sample=f"{p}_{g}", patient=p, condition=g, n_cells=len(ix)))
        return idxs, pd.DataFrame(meta).set_index("sample")
    X = a.layers["counts"].tocsr()
    allres, cand_rows, summ = {}, [], []
    ratio_p = pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv", sep="\t", index_col=0)
    ratio_g = pd.read_csv(f"{P}/results/04_source_filter/ratio_GSE164789.tsv", sep="\t", index_col=0)
    sff = pd.read_csv(f"{P}/results/04_source_filter/source_filter_final.tsv", sep="\t", index_col=0)
    for lin, labs in LINS.items():
        for variant in ["main", "aismia", "stageI"]:
            idxs, meta = units(labs, lin not in ("AT2L",), variant)
            if len(meta) < 4: log(f"[{lin}][{variant}] 配对不足"); continue
            counts = pd.DataFrame({k: np.asarray(X[v].sum(axis=0)).ravel() for k, v in idxs.items()}, index=a.var_names).T.astype(int)
            meta["condition"] = pd.Categorical(meta["condition"], categories=["Normal", "Precursor"])
            r = m14.deseq(counts, meta[["patient", "condition"]], "~ patient + condition", ["condition", "Precursor", "Normal"])
            r["p1_up"] = p1_of(r)
            lib = counts.div(counts.sum(axis=1), axis=0) * 1e6; lg = np.log1p(lib)
            amb = lg[[g for g in ["SFTPC", "SFTPB", "SFTPA1"] if g in lg.columns]].mean(axis=1)
            r["ambient_corr"] = lg.apply(lambda col: spearmanr(col, amb)[0] if col.std() > 0 else np.nan).reindex(r.index)
            allres[(lin, variant)] = r
            r.to_csv(f"{OUT}/deg_{lin.replace('/','_')}_{variant}.tsv", sep="\t"); meta.to_csv(f"{OUT}/units_{lin.replace('/','_')}_{variant}.tsv", sep="\t")
            if variant == "main": counts.to_csv(f"{OUT}/pseudobulk_counts_{lin.replace('/','_')}.tsv", sep="\t"); meta.to_csv(f"{OUT}/pseudobulk_meta_{lin.replace('/','_')}.tsv", sep="\t")
            sig = r[r["padj"] < 0.05]
            summ.append(dict(lineage=lin, variant=variant, n_pairs=len(meta) // 2, n_genes=int(r["padj"].notna().sum()), n_sig=len(sig),
                             up=int((sig.log2FoldChange > 0).sum()), down=int((sig.log2FoldChange < 0).sum())))
            log(f"[{lin}][{variant}] {summ[-1]}")
    pd.DataFrame(summ).to_csv(f"{OUT}/summary_deg.tsv", sep="\t", index=False)
    # ---------- 候选 ----------
    for lin in ["AT2L", "Myeloid", "Fibroblast", "T/NK", "Endothelial"]:
        r = allres[(lin, "main")].copy()
        r["is_prolif"] = r.index.isin(prolif2)
        if lin == "AT2L":
            r["flag_source"] = sff["final_flag"].reindex(r.index).fillna(False).astype(bool)
        else:
            L = {"Myeloid": "Myeloid", "Fibroblast": "Fibroblast", "T/NK": "T/NK", "Endothelial": "Endothelial"}[lin]
            r["ratio_lin_to_epi_peng"] = (ratio_p[f"med_{L}"] / ratio_p["epi"].replace(0, np.nan)).reindex(r.index)
            fl = r["ratio_lin_to_epi_peng"] < 0.2
            resc = fl & ((ratio_g[f"med_{L}"] / ratio_g["epi"].replace(0, np.nan)).reindex(r.index) >= 0.2)
            r["flag_source"] = (fl & ~resc) | r.index.isin(m14.AMBIENT) | (r["ambient_corr"] > 0.6)
        cand = r[(r["padj"] < 0.05) & (r["log2FoldChange"] > 0)].copy()
        cand["pass_source"] = ~cand["flag_source"]; cand["pass_nonprolif"] = ~cand["is_prolif"]
        cand["lineage"] = lin
        for v in ["aismia", "stageI"]:
            rv = allres.get((lin, v))
            cand[f"retained_{v}"] = False if rv is None else ((rv["padj"].reindex(cand.index) < 0.05) & (rv["log2FoldChange"].reindex(cand.index) > 0))
        cand["eligible"] = cand["pass_source"] & cand["pass_nonprolif"]
        cand.to_csv(f"{OUT}/candidates_{lin.replace('/','_')}_all_up.tsv", sep="\t")
        cand_rows.append(cand)
        e = cand[cand["eligible"]]
        log(f"[{lin}] 前体>正常 FDR<0.05 {len(cand)}；增殖基因 {int(cand.is_prolif.sum())}；来源过滤剔除 {int((~cand.pass_source).sum())}；合格候选 {len(e)}；其中 AIS/MIA 敏感性保留 {int(e.retained_aismia.sum())}，排除 II–III 期保留 {int(e.retained_stageI.sum())}")
    pd.concat(cand_rows).to_csv(f"{OUT}/candidates_stage2_all.tsv", sep="\t")
    pd.concat([c[c.eligible] for c in cand_rows]).to_csv(f"{OUT}/candidates_stage2_eligible.tsv", sep="\t")

    # ---------- 细胞水平：KAC、CAS（AT2L）----------
    log("计算细胞水平评分（本研究自定义阈值）")
    m_at2 = o["fl"].isin(LINS["AT2L"]).values
    s2 = a[m_at2].copy(); s2.obs = s2.obs[["patient", "stage"]].copy()
    for nm, genes in [("KAC", kac), ("CAS", CAS), ("KAC_noSPINK1", [g for g in kac if g != "SPINK1"])]:
        pres = [g for g in genes if g in s2.var_names]
        sc.tl.score_genes(s2, pres, score_name=nm, random_state=0)
        log(f"  {nm}：{len(pres)}/{len(genes)} 个基因可用")
    nrm_mask = (s2.obs["stage"] == "Normal").values
    thr = {nm: np.percentile(s2.obs.loc[nrm_mask, nm], 90) for nm in ["KAC", "CAS"]}
    log(f"  阈值（Normal AT2 谱系细胞第 90 百分位）：{thr}")
    idx_main, meta_main = units(LINS["AT2L"], False, "main")
    rows = []
    # 重新用 s2 的行索引：units 返回的是 a 的位置；转为 s2 的位置
    pos_map = -np.ones(a.n_obs, dtype=int); pos_map[np.where(m_at2)[0]] = np.arange(m_at2.sum())
    for k, ix in idx_main.items():
        sidx = pos_map[ix]
        d = dict(sample=k, patient=k.split("_")[0], condition=k.split("_")[1], n_cells=len(ix))
        d["KAC_frac"] = float((s2.obs["KAC"].values[sidx] > thr["KAC"]).mean()); d["CAS_frac"] = float((s2.obs["CAS"].values[sidx] > thr["CAS"]).mean())
        d["KAC_mean"] = float(s2.obs["KAC"].values[sidx].mean()); d["CAS_mean"] = float(s2.obs["CAS"].values[sidx].mean())
        rows.append(d)
    fr = pd.DataFrame(rows).set_index("sample")
    # 相关性
    cell_corr = {}
    for lab, mk in [("全部", np.ones(s2.n_obs, bool)), ("Normal", nrm_mask), ("前体", np.isin(s2.obs["stage"].values, ["AAH", "AIS", "MIA"]))]:
        cell_corr[lab] = (spearmanr(s2.obs.loc[mk, "CAS"], s2.obs.loc[mk, "KAC"])[0], spearmanr(s2.obs.loc[mk, "CAS"], s2.obs.loc[mk, "KAC_noSPINK1"])[0], int(mk.sum()))
    del s2
    # ---------- 巨噬细胞：IL1B+、SIRPA+ ----------
    m_mac = (o["fl"] == "Macrophage").values & (epcam < 2)
    mm = a[m_mac]; il1b = np.asarray(mm[:, "IL1B"].X.todense()).ravel() if sp.issparse(mm.X) else np.asarray(mm[:, "IL1B"].X).ravel()
    sirpa = np.asarray(mm[:, "SIRPA"].layers["counts"].todense()).ravel()
    st_m = mm.obs["stage"].astype(str).values; pt_m = mm.obs["patient"].astype(str).values
    thr_il1b = np.percentile(il1b[st_m == "Normal"], 90)
    il1b_pos = (il1b > thr_il1b) if thr_il1b > 0 else (il1b > 0)
    log(f"  IL1B 第 90 百分位（Normal 巨噬细胞）= {thr_il1b:.4f}；IL1B+ 定义：{'> 阈值' if thr_il1b>0 else 'IL1B>0'}；Normal 中 IL1B+ 比例 {il1b_pos[st_m=='Normal'].mean()*100:.1f}%")
    mrows = []
    for p in sorted(set(pt_m), key=lambda x: int(x[1:])):
        for g, stg in [("Normal", ["Normal"]), ("Precursor", ["AAH", "AIS", "MIA"]), ("PrecursorAISMIA", ["AIS", "MIA"])]:
            sel = (pt_m == p) & np.isin(st_m, stg)
            if sel.sum() >= MIN_CELLS: mrows.append(dict(patient=p, group=g, n=int(sel.sum()), IL1B_pos_frac=float(il1b_pos[sel].mean()), SIRPA_pos_frac=float((sirpa[sel] >= 1).mean())))
    mf = pd.DataFrame(mrows)
    fr.to_csv(f"{OUT}/cell_fractions_AT2L.tsv", sep="\t"); mf.to_csv(f"{OUT}/cell_fractions_macrophage.tsv", sep="\t", index=False)

    # ---------- 样本水平 ssGSEA（Peng 方法；gseapy ssGSEA 默认参数）----------
    cnt = pd.read_csv(f"{OUT}/pseudobulk_counts_AT2L.tsv", sep="\t", index_col=0)
    ld = np.log1p(cnt.div(cnt.sum(axis=1), axis=0) * 1e6).T
    gs = {"KAC": [g for g in kac if g in ld.index], "NFKB": [g for g in nfkb if g in ld.index]}
    ss = gp.ssgsea(data=ld, gene_sets=gs, outdir=None, no_plot=True, threads=4, min_size=15)
    sres = ss.res2d.pivot(index="Name", columns="Term", values="NES").astype(float)
    sres.to_csv(f"{OUT}/ssgsea_AT2L_pseudobulk.tsv", sep="\t")

    # ---------- 假说检验 ----------
    def paired_stat(fr_, col, group_a="Precursor", group_b="Normal", pat_col="patient", cond_col="condition", pairs_variant=None):
        a_ = fr_[fr_[cond_col] == group_a].set_index(pat_col)[col]; b_ = fr_[fr_[cond_col] == group_b].set_index(pat_col)[col]
        d = (a_ - b_.reindex(a_.index)).dropna()
        if pairs_variant is not None: d = d[~d.index.isin(pairs_variant)]
        return d
    tests = []
    def add(fam, tid, primary, desc, stat, n, p_up, p2=np.nan, extra=""):
        tests.append(dict(family=fam, test=tid, primary=primary, description=desc, statistic=stat, n=n, p_one_sided=p_up, p_two_sided=p2, note=extra))
    def deseq_row(lin, gene):
        r = allres[(lin, "main")]
        return (r.loc[gene, "log2FoldChange"], r.loc[gene, "pvalue"], r.loc[gene, "padj"], r.loc[gene, "p1_up"]) if gene in r.index and pd.notna(r.loc[gene, "pvalue"]) else (np.nan,) * 4
    n_main = len(allres[("AT2L", "main")]) and int(len(pd.read_csv(f"{OUT}/units_AT2L_main.tsv", sep="\t")) / 2)
    for tid, lin, gene, primary in [("H1a", "AT2L", "IL1R1", True), ("H1b", "Myeloid", "IL1B", False), ("H1c-RELA", "AT2L", "RELA", False), ("H1c-RELB", "AT2L", "RELB", False), ("H1c-NFKB1", "AT2L", "NFKB1", False)]:
        l, p2, q, p1 = deseq_row(lin, gene); add("H1", tid, primary, f"{lin} {gene} 伪批量表达（前体 vs Normal）", l, n_main, p1, p2, f"log2FC；双侧 FDR={q:.3g}")
    d = paired_stat(fr, "KAC_frac"); add("H1", "H1d", False, "KAC 细胞比例（AT2 谱系；本研究自定义阈值）", float(d.mean()), len(d), wilcoxon(d, alternative="greater").pvalue, wilcoxon(d).pvalue, f"均值差；前体更高 {int((d>0).sum())}/{len(d)}")
    dm = mf[mf.group == "Precursor"].set_index("patient")["IL1B_pos_frac"] - mf[mf.group == "Normal"].set_index("patient")["IL1B_pos_frac"].reindex(mf[mf.group == "Precursor"].set_index("patient").index)
    dm = dm.dropna(); add("H1", "H1e", True, "IL1B+ 巨噬细胞比例（本研究自定义阈值）", float(dm.mean()), len(dm), wilcoxon(dm, alternative="greater").pvalue, wilcoxon(dm).pvalue, f"均值差；前体更高 {int((dm>0).sum())}/{len(dm)}")
    for tid, lin, gene in [("H2a", "AT2L", "CD47"), ("H2b", "Macrophage", "SIRPA")]:
        l, p2, q, p1 = deseq_row(lin, gene); add("H2", tid, False, f"{lin} {gene} 伪批量表达", l, int(len(pd.read_csv(f'{OUT}/units_{lin}_main.tsv', sep='\t')) / 2), p1, p2, f"log2FC；双侧 FDR={q:.3g}")
    ds = mf[mf.group == "Precursor"].set_index("patient")["SIRPA_pos_frac"] - mf[mf.group == "Normal"].set_index("patient")["SIRPA_pos_frac"].reindex(mf[mf.group == "Precursor"].set_index("patient").index)
    ds = ds.dropna(); add("H2", "H2c", False, "SIRPA+ 巨噬细胞比例（SIRPA 计数>=1）", float(ds.mean()), len(ds), wilcoxon(ds, alternative="greater").pvalue, wilcoxon(ds).pvalue, f"均值差；前体更高 {int((ds>0).sum())}/{len(ds)}")
    dc_ = paired_stat(fr, "CAS_frac"); add("H3", "H3", False, "CAS 样上皮比例（AT2 谱系；本研究自定义评分与阈值）", float(dc_.mean()), len(dc_), wilcoxon(dc_, alternative="greater").pvalue, wilcoxon(dc_).pvalue, f"均值差；前体更高 {int((dc_>0).sum())}/{len(dc_)}")
    T = pd.DataFrame(tests)
    T["p_holm_family"] = np.nan
    for fam, g in T.groupby("family"):
        T.loc[g.index, "p_holm_family"] = multipletests(g["p_one_sided"], method="holm")[1]
    prim = T[T.primary]; T["p_holm_primary_only"] = np.nan
    T.loc[prim.index, "p_holm_primary_only"] = multipletests(prim["p_one_sided"], method="holm")[1]
    T["supported"] = (T["statistic"] > 0) & (T["p_one_sided"] < 0.05) & (T["p_holm_family"] < 0.05)
    T.to_csv(f"{OUT}/hypotheses_H1_H2_H3.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 60)
    log("假说检验：\n" + T[["family", "test", "primary", "description", "statistic", "n", "p_one_sided", "p_holm_family", "supported", "note"]].round(5).to_string(index=False))
    # ssGSEA 样本水平（只报告未校正 P）
    srows = []
    for term in sres.columns:
        x = sres[term]; pre = x[[i for i in x.index if i.endswith("_Precursor")]]; nm = x[[i for i in x.index if i.endswith("_Normal")]]
        pre.index = [i.split("_")[0] for i in pre.index]; nm.index = [i.split("_")[0] for i in nm.index]
        dd = (pre - nm.reindex(pre.index)).dropna()
        srows.append(dict(score=f"ssGSEA {term}", n=len(dd), mean_diff=dd.mean(), p_one_sided=wilcoxon(dd, alternative="greater").pvalue, up_in=int((dd > 0).sum())))
    pd.DataFrame(srows).to_csv(f"{OUT}/ssgsea_paired_tests.tsv", sep="\t", index=False)
    log("ssGSEA 样本水平（未校正单侧 P）：" + "; ".join(f"{r['score']}: 均值差 {r['mean_diff']:+.3f}，P={r['p_one_sided']:.4g}，前体更高 {r['up_in']}/{r['n']}" for r in srows))
    # H3 相关
    sc_ = fr.copy(); sc_["g"] = sc_["condition"]
    rho_s = spearmanr(sc_["CAS_frac"], sc_["KAC_frac"])[0]
    dd = pd.DataFrame({"dCAS": paired_stat(fr, "CAS_frac"), "dKAC": paired_stat(fr, "KAC_frac")}).dropna(); rho_d = spearmanr(dd["dCAS"], dd["dKAC"])[0]
    log("H3 相关（Spearman）：细胞水平 CAS×KAC（含 SPINK1 / 去掉 SPINK1 的 KAC；细胞数）：" + "; ".join(f"{k}: {v[0]:.3f} / {v[1]:.3f}（{v[2]}）" for k, v in cell_corr.items()))
    log(f"H3 相关：样本水平 CAS 比例×KAC 比例（Normal+前体伪批量，n={len(sc_)}）rho={rho_s:.3f}；前体−Normal 差值相关（n={len(dd)}）rho={rho_d:.3f}")
    pd.DataFrame([dict(level="cell_all", rho_with_SPINK1=cell_corr["全部"][0], rho_without_SPINK1=cell_corr["全部"][1], n=cell_corr["全部"][2]),
                  dict(level="cell_Normal", rho_with_SPINK1=cell_corr["Normal"][0], rho_without_SPINK1=cell_corr["Normal"][1], n=cell_corr["Normal"][2]),
                  dict(level="cell_Precursor", rho_with_SPINK1=cell_corr["前体"][0], rho_without_SPINK1=cell_corr["前体"][1], n=cell_corr["前体"][2]),
                  dict(level="sample_fraction", rho_with_SPINK1=rho_s, n=len(sc_)), dict(level="paired_difference", rho_with_SPINK1=rho_d, n=len(dd))]).to_csv(f"{OUT}/H3_correlations.tsv", sep="\t", index=False)
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
