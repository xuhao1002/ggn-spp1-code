#!/usr/bin/env python
"""
第 1 层外部验证（预设见 docs/项目进展与交接.md「外部验证的最终预设」）。
- 验证集：GSE164789（病灶 -T 样本）、GSE189357；检验 IAC vs (AIS+MIA)，pydeseq2，设计 ~ condition。
- 单位：病灶样本（患者×组织学）。GSE164789 同患者同组织学的多个 -T 样本合并（LRH-T + LRH-KT → LRH_MIA）；
  HYY 的 IAC 病灶（HYY-T1）与 MIA 病灶（HYY-T2）分别作为两个单位纳入；每个单位至少 50 个上皮细胞。
- 上皮细胞 = 注释为上皮（AT2、AT1、Club/secretory、Ciliated、Basal）的细胞；剔除 PTPRC 原始计数 >= 2 的上皮细胞（跨谱系共表达，与 01b 同规则）。
- 判定：参与校正的基因 = A 的 FDR<0.05 且通过来源过滤（含被救回）的基因；在每个验证集内，对「按预期方向的单侧 P」做 BH 校正；
  q<0.05 且方向一致才算通过；名义 P 同时报告；至少一个验证集通过 = 通过外部验证。
- 优先级分层（看到差异基因数后加的，不作入选条件）：梯队1 = |log2FC_A|>0.5 且 A∩B 同向显著；梯队2 = 满足其一；梯队3 = 其余。
- 敏感性（只报告，不参与判定）：GSE164789 用 ~ lesion_pair + TG（TG = A / T_前体 / T_IAC），检验 (T_IAC−A) − (T_前体−A)。
- 说明：验证集是患者间比较，通过表示方向稳健，不是严格复现；C 与 A 共享 16 对样本，不是独立验证；GSE189357 只有 3 vs 6，检验力很低。
追加（2026-10-08，事后发现不一致而加）：主要分析 A 用 AT2 谱系，而上面预设的验证集用了全部上皮。
  `python scripts/06_external_validation.py at2` = 验证集只取 AT2 + Club/secretory（对应 Peng 的 AT2 + AT2/Club-like；验证集无增殖上皮标签），输出 results/06_validation_at2lineage/。两套结果并列。
用法：nohup /opt/anaconda3/envs/ggn_target/bin/python scripts/06_external_validation.py [all|at2] > results/06_console.txt 2>&1 &
输出：results/06_validation/
"""
import os, sys, time, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

B = os.path.expanduser("~/全部数据库/GGN_target")
LINEAGE = sys.argv[1] if len(sys.argv) > 1 else "all"   # all=全部上皮（预设主分析）；at2=只取 AT2 谱系（事后追加的敏感性分析）
OUT = f"{B}/results/06_validation" + ("" if LINEAGE == "all" else "_at2lineage"); os.makedirs(OUT, exist_ok=True)
EPI = ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal"] if LINEAGE == "all" else ["AT2", "Club/secretory"]
MIN_CELLS = 50
LOG = open(f"{OUT}/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

def load(ds):
    d = f"{B}/results/05_validation/{ds}"
    a = sc.read_h5ad(f"{d}/{ds}_annotated.h5ad")
    ann = pd.read_csv(f"{d}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
    a.obs["final_label"] = a.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
    is_epi = a.obs["final_label"].isin(EPI).values
    ptprc = a[:, "PTPRC"].layers["counts"]
    ptprc = np.asarray(ptprc.todense()).ravel() if sp.issparse(ptprc) else np.asarray(ptprc).ravel()
    drop = is_epi & (ptprc >= 2)
    log(f"[{ds}] 上皮 {is_epi.sum()}；其中 PTPRC>=2 剔除 {int(drop.sum())}（{drop.sum()/is_epi.sum()*100:.1f}%）")
    keep = is_epi & ~drop
    return a[keep].copy()

def pseudobulk(a, unit_col, meta_cols, name):
    X = a.layers["counts"].tocsr()
    rows, labels, meta, excl = [], [], [], []
    for u, idx in a.obs.groupby(unit_col, observed=True).indices.items():
        if len(idx) < MIN_CELLS:
            excl.append(f"{u}({len(idx)})"); continue
        rows.append(np.asarray(X[idx].sum(axis=0)).ravel()); labels.append(u)
        m = a.obs.iloc[idx][meta_cols].iloc[0].to_dict(); m["n_cells"] = len(idx); meta.append(m)
    log(f"[{name}] 伪批量单位 {len(labels)} 个；<50 个上皮细胞被排除：{excl or '无'}")
    return pd.DataFrame(np.vstack(rows), index=labels, columns=a.var_names).astype(int), pd.DataFrame(meta, index=labels)

def filt(counts):
    n = len(counts)
    keep = ((counts >= 10).sum(axis=0) >= np.ceil(n / 3)) & (counts.sum(axis=0) >= 50)
    return counts.loc[:, keep.values]

def deseq_cond(counts, meta, name):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    c = filt(counts)
    m = meta[["condition"]].copy(); m["condition"] = pd.Categorical(m["condition"], categories=["Precursor", "IAC"])
    dds = DeseqDataSet(counts=c, metadata=m, design="~ condition", refit_cooks=True, n_cpus=6, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "IAC", "Precursor"], quiet=True); ds.summary()
    r = ds.results_df.copy(); r.index.name = "gene"
    r.to_csv(f"{OUT}/deg_{name}.tsv", sep="\t")
    log(f"[{name}] 检测基因 {c.shape[1]}；单位 {len(c)}（IAC {int((m['condition']=='IAC').sum())}，前体 {int((m['condition']=='Precursor').sum())}）；"
        f"名义 P<0.05 的基因 {int((r['pvalue']<0.05).sum())}")
    return r

def main():
    log("===== 开始 =====")
    A = pd.read_csv(f"{B}/results/03_pseudobulk/deg_A_main.tsv", sep="\t", index_col=0)
    Bst = pd.read_csv(f"{B}/results/03_pseudobulk/deg_B_state.tsv", sep="\t", index_col=0)
    S = pd.read_csv(f"{B}/results/04_source_filter/source_filter_final.tsv", sep="\t", index_col=0)
    Asig = A[A["padj"] < 0.05].copy()
    Asig = Asig.join(S[["flag_peng", "rescued", "final_flag", "lineage_peng", "ratio_peng", "ratio_GSE164789"]], how="left")
    Asig["final_flag"] = Asig["final_flag"].fillna(False).astype(bool)
    Asig["pass_source"] = ~Asig["final_flag"]
    Asig["B_lfc"] = Bst["log2FoldChange"].reindex(Asig.index); Asig["B_padj"] = Bst["padj"].reindex(Asig.index)
    Asig["B_same_sig"] = (Asig["B_padj"] < 0.05) & (np.sign(Asig["B_lfc"]) == np.sign(Asig["log2FoldChange"]))
    c1, c2 = Asig["log2FoldChange"].abs() > 0.5, Asig["B_same_sig"]
    Asig["tier"] = np.where(c1 & c2, 1, np.where(c1 | c2, 2, 3))
    cand = Asig[Asig["pass_source"]]
    log(f"A 显著基因 {len(Asig)}；通过来源过滤 {len(cand)}")

    results = {}
    # ---- GSE164789 ----
    a = load("GSE164789")
    o = a.obs
    o["hist"] = o["histology"].astype(str)
    t = a[(o["tissue"] == "T").values].copy()
    t.obs["unit"] = t.obs["patient"].astype(str) + "_" + t.obs["histology"].astype(str)
    t.obs["condition"] = np.where(t.obs["histology"].astype(str) == "IAC", "IAC", "Precursor")
    cnt, meta = pseudobulk(t, "unit", ["condition", "patient", "histology"], "GSE164789")
    meta.to_csv(f"{OUT}/units_GSE164789.tsv", sep="\t")
    results["GSE164789"] = deseq_cond(cnt, meta, "GSE164789")

    # 敏感性：病灶 vs 癌旁的交互（只报告）
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        a.obs["pair"] = a.obs["patient"].astype(str) + "_" + a.obs["histology"].astype(str)
        a.obs["TG"] = np.where(a.obs["tissue"] == "A", "A", np.where(a.obs["histology"].astype(str) == "IAC", "T_IAC", "T_Pre"))
        a.obs["unit2"] = a.obs["pair"] + "|" + a.obs["TG"]
        cnt2, meta2 = pseudobulk(a, "unit2", ["pair", "TG"], "GSE164789_interaction")
        ok = meta2.groupby("pair").size() == 2          # 病灶与癌旁都 >=50 的配对
        meta2 = meta2[meta2["pair"].isin(ok.index[ok])]; cnt2 = cnt2.loc[meta2.index]
        c = filt(cnt2)
        m = meta2[["pair", "TG"]].copy(); m["TG"] = pd.Categorical(m["TG"], categories=["A", "T_Pre", "T_IAC"])
        dds = DeseqDataSet(counts=c, metadata=m, design="~ pair + TG", refit_cooks=True, n_cpus=6, quiet=True)
        dds.deseq2()
        cols = list(dds.obsm["design_matrix"].columns)
        v = np.zeros(len(cols)); v[cols.index("TG[T.T_IAC]")] = 1; v[cols.index("TG[T.T_Pre]")] = -1
        ds = DeseqStats(dds, contrast=v, quiet=True); ds.summary()
        r = ds.results_df.copy(); r.index.name = "gene"; r.to_csv(f"{OUT}/deg_GSE164789_interaction.tsv", sep="\t")
        results["interaction"] = r
        log(f"[交互项] 配对 {len(m)//2} 个；检测基因 {c.shape[1]}；名义 P<0.05 {int((r['pvalue']<0.05).sum())}；FDR<0.05 {int((r['padj']<0.05).sum())}")
    except Exception as e:
        log(f"[交互项] 失败（不影响判定）：{e!r}")
    del a, t

    # ---- GSE189357 ----
    b = load("GSE189357")
    b.obs["condition"] = np.where(b.obs["histology"].astype(str) == "IAC", "IAC", "Precursor")
    cnt, meta = pseudobulk(b, "sample", ["condition", "patient", "histology"], "GSE189357")
    meta.to_csv(f"{OUT}/units_GSE189357.tsv", sep="\t")
    results["GSE189357"] = deseq_cond(cnt, meta, "GSE189357")

    # ---- 判定 ----
    out = cand.copy()
    for ds in ["GSE164789", "GSE189357"]:
        r = results[ds].reindex(out.index)
        out[f"{ds}_lfc"] = r["log2FoldChange"]; out[f"{ds}_p"] = r["pvalue"]
        exp = np.sign(out["log2FoldChange"])
        same = np.sign(r["log2FoldChange"]) == exp
        p1 = np.where(same, r["pvalue"] / 2, 1 - r["pvalue"] / 2)
        out[f"{ds}_same_dir"] = same
        out[f"{ds}_p1"] = np.where(r["pvalue"].notna(), p1, np.nan)
        ok = out[f"{ds}_p1"].notna()
        q = np.full(len(out), np.nan)
        q[ok.values] = multipletests(out.loc[ok, f"{ds}_p1"], method="fdr_bh")[1]
        out[f"{ds}_q"] = q
        out[f"{ds}_pass"] = same & (out[f"{ds}_q"] < 0.05)
        log(f"[{ds}] 参与校正基因 {int(ok.sum())}（A∩来源过滤共 {len(out)}）；方向一致 {int((same & ok).sum())}（{(same & ok).sum()/ok.sum()*100:.1f}%）；"
            f"名义单侧 P<0.05 且同向 {int((same & (out[f'{ds}_p1']<0.05)).sum())}；q<0.05 通过 {int(out[f'{ds}_pass'].sum())}")
    out["pass_any"] = out["GSE164789_pass"] | out["GSE189357_pass"]
    out["pass_both"] = out["GSE164789_pass"] & out["GSE189357_pass"]
    if "interaction" in results:
        out["interaction_lfc"] = results["interaction"]["log2FoldChange"].reindex(out.index)
        out["interaction_padj"] = results["interaction"]["padj"].reindex(out.index)
    out = out.sort_values(["pass_any", "tier", "padj"], ascending=[False, True, True])
    out.to_csv(f"{OUT}/validation_results.tsv", sep="\t")
    cp = out[out["pass_any"]]
    cp.to_csv(f"{OUT}/layer1_candidates.tsv", sep="\t")
    log(f"第 1 层通过外部验证（至少一个验证集 q<0.05 且同向）：{len(cp)}；两个都通过 {int(out['pass_both'].sum())}")
    log("梯队分布（通过者）：" + str(cp["tier"].value_counts().sort_index().to_dict()) + "；（全部候选）" + str(out["tier"].value_counts().sort_index().to_dict()))
    top = cp.sort_values(["tier", "padj"]).head(30)
    log("通过者前 30（按梯队、A 的 padj）：" + ", ".join(f"{g}(T{int(x.tier)},{x.log2FoldChange:+.1f})" for g, x in top.iterrows()))
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
