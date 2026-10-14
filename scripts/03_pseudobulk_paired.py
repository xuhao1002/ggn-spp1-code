#!/usr/bin/env python
"""
第三步：同一患者内「前体病变 vs LUAD」的配对伪批量差异分析（pydeseq2，设计 ~ patient + condition）。

用法：
    cd ~/全部数据库/GGN_target && nohup /opt/anaconda3/envs/ggn_target/bin/python scripts/03_pseudobulk_paired.py > results/03_console.txt 2>&1 &

预设（2026-10-08 运行前写进 docs/项目进展与交接.md）：
- 伪批量 = 患者 × (Precursor | LUAD)；前体 = 同患者所有 AAH/AIS/MIA 样本合并；LUAD = 同患者所有 LUAD 样本合并；原始计数相加。
- 对比：LUAD vs Precursor；FDR<0.05（pydeseq2 的 BH 校正）；每个伪批量至少 50 个细胞，否则该患者排除。
- 基因过滤：原始计数 >= 10 的伪批量样本数 >= 总样本数的 1/3，且总计数 >= 50。
- 分析：
    A_main        AT2 谱系（AT2 + AT2/Club-like epithelial + Proliferating epithelial），全部配对（22 对）
    B_state       仅 AT2/Club-like epithelial 亚类（区分细胞状态转变 vs 状态内表达改变）
    C_AISMIA      A 的敏感性：前体仅 AIS/MIA（P9、P22 只保留其 AIS 细胞）
    D_allepi      A 的敏感性：全部上皮（final_compartment == Epithelial）
- 汇总：各分析上/下调基因数；A 的前 50 个基因；A∩B 同向显著基因；C 与 A 的方向一致率（以 A 中 FDR<0.05 基因为准）。
输出：results/03_pseudobulk/
    pseudobulk_counts_{analysis}.tsv、deg_{analysis}.tsv、deg_A_top50.tsv、overlap_A_B.tsv、summary.tsv、run_log.txt
"""
import os, json, time, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
warnings.filterwarnings("ignore")

BASE = os.path.expanduser("~/全部数据库/GGN_target")
OUT = os.path.join(BASE, "results/03_pseudobulk")
os.makedirs(OUT, exist_ok=True)
LIN = ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"]
MIN_CELLS = 50
LOG = open(os.path.join(OUT, "run_log.txt"), "a", encoding="utf-8")
T0 = time.time()
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

def pseudobulk(a, mask, name):
    """返回 (counts[样本×基因], metadata, 排除说明)"""
    o = a.obs[["patient", "stage", "sample"]].copy()
    o["grp"] = o["stage"].map(lambda s: "LUAD" if s == "LUAD" else ("Precursor" if s in ("AAH", "AIS", "MIA") else "Normal"))
    keep = mask & (o["grp"] != "Normal").values
    o = o[keep]
    X = a.layers["counts"][np.where(keep)[0]]
    rows, labels, meta, excl = [], [], [], []
    for p in sorted(o["patient"].unique(), key=lambda x: int(x[1:])):
        cells = {g: np.where(((o["patient"] == p) & (o["grp"] == g)).values)[0] for g in ["Precursor", "LUAD"]}
        n = {g: len(v) for g, v in cells.items()}
        if min(n.values()) < MIN_CELLS:
            excl.append(f"{p}(前体 {n['Precursor']}，LUAD {n['LUAD']})"); continue
        for g in ["Precursor", "LUAD"]:
            rows.append(np.asarray(X[cells[g]].sum(axis=0)).ravel())
            labels.append(f"{p}_{g}")
            meta.append(dict(patient=p, condition=g, n_cells=n[g]))
    counts = pd.DataFrame(np.vstack(rows), index=labels, columns=a.var_names).astype(int)
    meta = pd.DataFrame(meta, index=labels)
    log(f"[{name}] 伪批量：{len(meta)//2} 对；排除：{excl or '无'}")
    return counts, meta, excl

def run_deseq(counts, meta, name):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    n = len(counts)
    keep = ((counts >= 10).sum(axis=0) >= np.ceil(n / 3)) & (counts.sum(axis=0) >= 50)
    c = counts.loc[:, keep.values]
    log(f"[{name}] 基因过滤后 {c.shape[1]} / {counts.shape[1]}")
    m = meta[["patient", "condition"]].copy()
    m["condition"] = pd.Categorical(m["condition"], categories=["Precursor", "LUAD"])
    dds = DeseqDataSet(counts=c, metadata=m, design="~ patient + condition", refit_cooks=True, n_cpus=6, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "LUAD", "Precursor"], quiet=True)
    ds.summary()
    r = ds.results_df.copy()
    r.index.name = "gene"
    r = r.sort_values("padj")
    r.to_csv(os.path.join(OUT, f"deg_{name}.tsv"), sep="\t")
    counts.to_csv(os.path.join(OUT, f"pseudobulk_counts_{name}.tsv"), sep="\t")
    meta.to_csv(os.path.join(OUT, f"pseudobulk_meta_{name}.tsv"), sep="\t")
    return r

def main():
    log("===== 开始 =====")
    a = sc.read_h5ad(os.path.join(BASE, "results/01_snRNA/epithelial_cells.h5ad"))
    log(f"上皮细胞：{a.shape}")
    fl = a.obs["final_label"].astype(str).values
    st = a.obs["stage"].astype(str).values
    analyses = {
        "A_main": np.isin(fl, LIN),
        "B_state": fl == "AT2/Club-like epithelial",
        "C_AISMIA": np.isin(fl, LIN) & ~(st == "AAH"),
        "D_allepi": np.ones(a.n_obs, dtype=bool),
    }
    res, summ = {}, []
    for name, mask in analyses.items():
        counts, meta, excl = pseudobulk(a, mask, name)
        r = run_deseq(counts, meta, name)
        sig = r[r["padj"] < 0.05]
        s = dict(analysis=name, n_pairs=len(meta) // 2, n_genes_tested=int(r["padj"].notna().sum()),
                 n_sig=len(sig), n_up_in_LUAD=int((sig["log2FoldChange"] > 0).sum()),
                 n_down_in_LUAD=int((sig["log2FoldChange"] < 0).sum()),
                 n_sig_abslfc_gt_0p5=int((sig["log2FoldChange"].abs() > 0.5).sum()), excluded=";".join(excl))
        summ.append(s); res[name] = r
        log(f"[{name}] 显著 {s['n_sig']}（上调 {s['n_up_in_LUAD']} / 下调 {s['n_down_in_LUAD']}）；|log2FC|>0.5 的 {s['n_sig_abslfc_gt_0p5']}")
    pd.DataFrame(summ).to_csv(os.path.join(OUT, "summary.tsv"), sep="\t", index=False)

    A = res["A_main"]
    A.head(50).to_csv(os.path.join(OUT, "deg_A_top50.tsv"), sep="\t")

    # A∩B
    B = res["B_state"]
    sa, sb = A[A["padj"] < 0.05], B[B["padj"] < 0.05]
    both = sa.join(sb[["log2FoldChange", "padj"]], rsuffix="_B", how="inner")
    both["same_direction"] = np.sign(both["log2FoldChange"]) == np.sign(both["log2FoldChange_B"])
    both.to_csv(os.path.join(OUT, "overlap_A_B.tsv"), sep="\t")
    n_same = int(both["same_direction"].sum())
    log(f"A∩B：A 显著 {len(sa)}，B 显著 {len(sb)}，两者都显著 {len(both)}，其中方向一致 {n_same}")

    # 一致率
    cons = []
    for nm in ["B_state", "C_AISMIA", "D_allepi"]:
        R = res[nm]
        j = sa.join(R[["log2FoldChange", "padj"]], rsuffix="_other", how="left")
        ok = j["log2FoldChange_other"].notna()
        same = (np.sign(j.loc[ok, "log2FoldChange"]) == np.sign(j.loc[ok, "log2FoldChange_other"]))
        common = A.index.intersection(R.index)
        rho = A.loc[common, "log2FoldChange"].corr(R.loc[common, "log2FoldChange"], method="spearman")
        both_sig = int((j.loc[ok, "padj_other"] < 0.05).sum())
        cons.append(dict(vs=nm, n_A_sig=len(sa), n_A_sig_tested_in_other=int(ok.sum()),
                         direction_concordance=round(float(same.mean()), 4), n_also_sig=both_sig,
                         spearman_lfc_all_common_genes=round(float(rho), 4)))
        log(f"A vs {nm}：A 显著基因中方向一致率 {same.mean():.3f}（{int(same.sum())}/{int(ok.sum())}）；也显著 {both_sig}；全体共同基因 log2FC Spearman {rho:.3f}")
    pd.DataFrame(cons).to_csv(os.path.join(OUT, "concordance_vs_A.tsv"), sep="\t", index=False)
    log("A 前 20：" + ", ".join(f"{g}({r['log2FoldChange']:+.2f})" for g, r in A.head(20).iterrows()))
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
