#!/usr/bin/env python
"""
微环境第 1 层（标准见 docs/项目进展与交接.md「补充分析判定标准」，2026-10-08 用户确认，运行前写定）。
- 细胞群：T/NK、Myeloid、B/Plasma、Fibroblast、Endothelial（血管）、SMC/pericyte；每群各自做 22 对配对伪批量（~ 患者 + 病变类型，LUAD vs 前体，
  pydeseq2，FDR<0.05；每侧>=50 个细胞，基因过滤同上皮：原始计数>=10 的样本数>=1/3 且总计数>=50）。
- 所有非上皮细胞再剔除 EPCAM 原始计数>=2 的细胞（跨谱系共表达，与 01b 同规则）。
- 反向来源过滤：(a) 群内中位表达/上皮中位表达 <0.2 → 疑似上皮来源（用 04 的「谱系×样本」中位表达表；GSE164789 同法核对，Peng 标记而 GSE164789 比值>=0.2 的救回）；
  (b) 固定环境 RNA 名单；(c) 与伪批量内「SFTPC+SFTPB+SFTPA1 平均表达」的 Spearman 相关>0.6 → 环境 RNA 相关。
- 外部验证：GSE164789（病灶单位，IAC vs AIS+MIA）和 GSE189357，同一细胞群，pydeseq2 ~ condition，单侧 P 做 BH（校正基因=该群显著且通过来源过滤的基因），q<0.05 且同向；至少一个验证集通过。
输出 results/14_microenv/
"""
import os, sys, time, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/14_microenv"; os.makedirs(OUT, exist_ok=True)
LIN = {"T/NK": ["T cell", "NK", "Proliferating T/NK"],
       "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"],
       "B/Plasma": ["B cell", "Plasma", "Proliferating B/Plasma"],
       "Fibroblast": ["Fibroblast"], "Endothelial": ["Endothelial"], "SMC/pericyte": ["SMC/pericyte"]}
AMBIENT = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "NAPSA", "SCGB1A1", "SCGB3A1", "SCGB3A2", "SLPI", "WFDC2", "MUC1", "LAMP3", "AGER"]
MIN_CELLS = 50
LOG = open(f"{OUT}/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

def filt(c):
    n = len(c); k = ((c >= 10).sum(axis=0) >= np.ceil(n / 3)) & (c.sum(axis=0) >= 50); return c.loc[:, k.values]

def deseq(counts, meta, design, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    c = filt(counts)
    dds = DeseqDataSet(counts=c, metadata=meta, design=design, refit_cooks=True, n_cpus=6, quiet=True); dds.deseq2()
    ds = DeseqStats(dds, contrast=contrast, quiet=True); ds.summary()
    r = ds.results_df.copy(); r.index.name = "gene"; return r

def sums(a, idx_groups):
    X = a.layers["counts"].tocsr(); out = {}
    for k, idx in idx_groups.items(): out[k] = np.asarray(X[idx].sum(axis=0)).ravel()
    return out

def main():
    log("===== 开始 =====")
    ratio_p = pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv", sep="\t", index_col=0)
    ratio_g = pd.read_csv(f"{P}/results/04_source_filter/ratio_GSE164789.tsv", sep="\t", index_col=0)
    # ---------- Peng 伪批量与差异分析 ----------
    a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
    a.obs["lin"] = ""
    for L, labs in LIN.items(): a.obs.loc[a.obs["final_label"].astype(str).isin(labs), "lin"] = L
    epcam = np.asarray(a[:, "EPCAM"].layers["counts"].todense()).ravel()
    keep = (a.obs["lin"] != "") & (epcam < 2)
    log(f"非上皮六群细胞 {int((a.obs['lin']!='').sum())}；其中 EPCAM>=2 剔除 {int(((a.obs['lin']!='') & (epcam>=2)).sum())}")
    a = a[keep.values].copy()
    a.obs["grp"] = a.obs["stage"].astype(str).map(lambda s: "LUAD" if s == "LUAD" else ("Precursor" if s in ("AAH", "AIS", "MIA") else "Normal"))
    results = {}
    ambient_present = [g for g in ["SFTPC", "SFTPB", "SFTPA1"] if g in a.var_names]
    summ = []
    for L in LIN:
        sub = a[(a.obs["lin"] == L).values & (a.obs["grp"] != "Normal").values]
        groups, meta, excl = {}, [], []
        for p in sorted(sub.obs["patient"].unique(), key=lambda x: int(x[1:])):
            ok = True; idxs = {}
            for g in ["Precursor", "LUAD"]:
                idx = np.where(((sub.obs["patient"] == p) & (sub.obs["grp"] == g)).values)[0]
                if len(idx) < MIN_CELLS: ok = False
                idxs[g] = idx
            if not ok: excl.append(p); continue
            for g in ["Precursor", "LUAD"]:
                groups[f"{p}_{g}"] = idxs[g]; meta.append(dict(sample=f"{p}_{g}", patient=p, condition=g, n_cells=len(idxs[g])))
        S = sums(sub, groups)
        counts = pd.DataFrame(S, index=sub.var_names).T.astype(int)
        meta = pd.DataFrame(meta).set_index("sample")
        meta["condition"] = pd.Categorical(meta["condition"], categories=["Precursor", "LUAD"])
        r = deseq(counts, meta[["patient", "condition"]], "~ patient + condition", ["condition", "LUAD", "Precursor"])
        # 环境 RNA 相关
        lib = counts.div(counts.sum(axis=1), axis=0) * 1e6; lg = np.log1p(lib)
        amb = lg[ambient_present].mean(axis=1)
        cor = lg.apply(lambda col: spearmanr(col, amb)[0] if col.std() > 0 else np.nan)
        r["ambient_corr"] = cor.reindex(r.index)
        r.to_csv(f"{OUT}/deg_{L.replace('/','_')}.tsv", sep="\t"); meta.to_csv(f"{OUT}/units_{L.replace('/','_')}.tsv", sep="\t")
        results[L] = r
        sig = r[r["padj"] < 0.05]
        summ.append(dict(lineage=L, n_pairs=len(meta) // 2, excluded=";".join(excl), n_genes=int(r["padj"].notna().sum()), n_sig=len(sig),
                         up=int((sig.log2FoldChange > 0).sum()), down=int((sig.log2FoldChange < 0).sum())))
        log(f"[{L}] {summ[-1]}")
    del a
    # ---------- 反向来源过滤 ----------
    flag_rows = []
    for L, r in results.items():
        r = r.copy()
        medL = ratio_p[f"med_{L}"] if f"med_{L}" in ratio_p else None
        r["ratio_lin_to_epi_peng"] = (medL / ratio_p["epi"].replace(0, np.nan)).reindex(r.index)
        r["flag_epi_origin_peng"] = r["ratio_lin_to_epi_peng"] < 0.2
        medG = ratio_g[f"med_{L}"] if f"med_{L}" in ratio_g else None
        r["ratio_lin_to_epi_gse164789"] = (medG / ratio_g["epi"].replace(0, np.nan)).reindex(r.index)
        r["rescued"] = r["flag_epi_origin_peng"] & (r["ratio_lin_to_epi_gse164789"] >= 0.2)
        r["flag_epi_origin"] = r["flag_epi_origin_peng"] & ~r["rescued"]
        r["flag_ambient_list"] = r.index.isin(AMBIENT)
        r["flag_ambient_corr"] = r["ambient_corr"] > 0.6
        r["pass_source"] = ~(r["flag_epi_origin"] | r["flag_ambient_list"] | r["flag_ambient_corr"])
        results[L] = r
        sig = r[r["padj"] < 0.05]
        log(f"[{L}] 显著 {len(sig)}；疑似上皮来源 {int(sig.flag_epi_origin.sum())}（救回 {int(sig.rescued.sum())}）；环境名单 {int(sig.flag_ambient_list.sum())}；环境相关 {int(sig.flag_ambient_corr.sum())}；通过来源过滤 {int(sig.pass_source.sum())}")
        r.to_csv(f"{OUT}/deg_{L.replace('/','_')}_flagged.tsv", sep="\t")
    # ---------- 外部验证 ----------
    val = {}
    for ds in ["GSE164789", "GSE189357"]:
        d = f"{P}/results/05_validation/{ds}"
        v = sc.read_h5ad(f"{d}/{ds}_annotated.h5ad")
        ann = pd.read_csv(f"{d}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
        v.obs["fl"] = v.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
        ep = np.asarray(v[:, "EPCAM"].layers["counts"].todense()).ravel()
        v.obs["lin"] = ""
        for L, labs in LIN.items(): v.obs.loc[v.obs["fl"].isin(labs), "lin"] = L
        v = v[((v.obs["lin"] != "") & (ep < 2)).values].copy()
        if ds == "GSE164789": v = v[(v.obs["tissue"] == "T").values].copy()
        v.obs["unit"] = v.obs["patient"].astype(str) + "_" + v.obs["histology"].astype(str) if ds == "GSE164789" else v.obs["sample"].astype(str)
        v.obs["condition"] = np.where(v.obs["histology"].astype(str) == "IAC", "IAC", "Precursor")
        for L in LIN:
            sub = v[(v.obs["lin"] == L).values]
            groups, meta = {}, []
            for u, idx in sub.obs.groupby("unit", observed=True).indices.items():
                if len(idx) < MIN_CELLS: continue
                groups[u] = idx; meta.append(dict(unit=u, condition=sub.obs["condition"].iloc[idx[0]], n_cells=len(idx)))
            if not meta:
                log(f"[{ds}][{L}] 没有单位达到 {MIN_CELLS} 个细胞，跳过"); continue
            meta = pd.DataFrame(meta).set_index("unit")
            if meta["condition"].nunique() < 2 or meta["condition"].value_counts().min() < 2:
                log(f"[{ds}][{L}] 单位不足，跳过 {meta['condition'].value_counts().to_dict()}"); continue
            counts = pd.DataFrame(sums(sub, groups), index=sub.var_names).T.astype(int)
            meta["condition"] = pd.Categorical(meta["condition"], categories=["Precursor", "IAC"])
            r = deseq(counts, meta[["condition"]], "~ condition", ["condition", "IAC", "Precursor"])
            val[(ds, L)] = r
            log(f"[{ds}][{L}] 单位 {meta['condition'].value_counts().to_dict()}；检测基因 {int(r['pvalue'].notna().sum())}")
        del v
    # ---------- 判定 ----------
    allrows = []
    for L, r in results.items():
        cand = r[(r["padj"] < 0.05) & r["pass_source"]].copy()
        cand["lineage"] = L
        for ds in ["GSE164789", "GSE189357"]:
            vr = val.get((ds, L))
            if vr is None: cand[f"{ds}_pass"] = False; continue
            vv = vr.reindex(cand.index)
            same = np.sign(vv["log2FoldChange"]) == np.sign(cand["log2FoldChange"])
            p1 = np.where(same, vv["pvalue"] / 2, 1 - vv["pvalue"] / 2)
            ok = vv["pvalue"].notna().values
            q = np.full(len(cand), np.nan)
            if ok.sum() > 0: q[ok] = multipletests(p1[ok], method="fdr_bh")[1]
            cand[f"{ds}_lfc"] = vv["log2FoldChange"]; cand[f"{ds}_p1"] = np.where(ok, p1, np.nan); cand[f"{ds}_q"] = q
            cand[f"{ds}_same_dir"] = same
            cand[f"{ds}_pass"] = same & (q < 0.05)
            log(f"[{L}][{ds}] 参与校正 {int(ok.sum())}/{len(cand)}；方向一致 {int((same & ok).sum())}（{(same&ok).sum()/max(ok.sum(),1)*100:.1f}%）；q<0.05 通过 {int(cand[f'{ds}_pass'].sum())}")
        cand["L1_pass"] = cand["GSE164789_pass"] | cand["GSE189357_pass"]
        allrows.append(cand)
        log(f"[{L}] 通过外部验证 {int(cand['L1_pass'].sum())}（上调 {int((cand['L1_pass'] & (cand.log2FoldChange>0)).sum())} / 下调 {int((cand['L1_pass'] & (cand.log2FoldChange<0)).sum())}）")
    al = pd.concat(allrows)
    al.to_csv(f"{OUT}/microenv_layer1_candidates.tsv", sep="\t")
    pd.DataFrame(summ).to_csv(f"{OUT}/summary_deg.tsv", sep="\t", index=False)
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
