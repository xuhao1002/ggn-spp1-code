#!/usr/bin/env python
"""
配体–受体分析（标准见 docs/项目进展与交接.md「补充分析判定标准」，2026-10-08 用户确认）。工具：LIANA+ 1.10.0（liana.mt.rank_aggregate，consensus 资源）。
- 候选相互作用：上皮 → 微环境群：配体 ∈ A 经验证的上调基因（预设 352 内上调 224），受体 ∈ 该群通过第 1 层外部验证的上调基因；
  微环境群 → 上皮：配体 ∈ 该群通过第 1 层外部验证的上调基因，受体 ∈ A 经验证的上调基因（复合体：任一亚基属于该集合）。
- 每个样本分别运行（每个细胞群>=20 个细胞，每群最多抽 500 个，随机种子 0），得分 = 1 − magnitude_rank（未检出=0；细胞群不足=缺失）。
- Peng（发现）：每个患者的前体样本（AAH/AIS/MIA）取平均、LUAD 样本取平均，配对 Wilcoxon（>=8 对），BH 校正；显著且 LUAD 更高 = 上升的相互作用。
- GSE164789（验证）：只用 -T 样本，同患者同组织学合并取平均，IAC vs (AIS+MIA) 单侧 Mann–Whitney，BH 校正（在 Peng 显著的相互作用中），q<0.05 且同向。
输出 results/17_ligand_receptor/
"""
import os, sys, time, warnings
import numpy as np, pandas as pd, scanpy as sc, liana as li
from joblib import Parallel, delayed
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/17_ligand_receptor"; os.makedirs(OUT, exist_ok=True)
LIN = {"T/NK": ["T cell", "NK", "Proliferating T/NK"], "Myeloid": ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"],
       "B/Plasma": ["B cell", "Plasma", "Proliferating B/Plasma"], "Fibroblast": ["Fibroblast"], "Endothelial": ["Endothelial"], "SMC/pericyte": ["SMC/pericyte"]}
EPI = ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal", "Proliferating epithelial", "AT2/Club-like epithelial", "AT2-like (patient-specific)"]
LOG = open(f"{OUT}/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

cand = pd.read_csv(f"{P}/results/07_candidates/candidates_683.tsv", sep="\t", index_col=0)
epi_up = set(cand.index[cand.in_main352 & (cand.direction == "up")])
mic = pd.read_csv(f"{P}/results/14_microenv/microenv_layer1_candidates.tsv", sep="\t", index_col=0)
mic = mic[mic["L1_pass"] & (mic["log2FoldChange"] > 0)]
mic_up = {L: set(mic.index[mic.lineage == L]) for L in LIN}
res = li.resource.select_resource("consensus")
res["lig"] = res["ligand"].str.split("_"); res["rec"] = res["receptor"].str.split("_")
pairs = []
for _, r in res.iterrows():
    for L in LIN:
        if any(g in epi_up for g in r["lig"]) and any(g in mic_up[L] for g in r["rec"]):
            pairs.append((r["ligand"], r["receptor"], "Epithelial", L))
        if any(g in mic_up[L] for g in r["lig"]) and any(g in epi_up for g in r["rec"]):
            pairs.append((r["ligand"], r["receptor"], L, "Epithelial"))
pairs = pd.DataFrame(pairs, columns=["ligand", "receptor", "source", "target"]).drop_duplicates()
sub_res = pairs[["ligand", "receptor"]].drop_duplicates()
log(f"候选相互作用：{len(pairs)}（含方向与细胞群），资源子集 {len(sub_res)} 对；按群：{pairs.groupby(['source','target']).size().to_dict()}")
pairs.to_csv(f"{OUT}/candidate_pairs.tsv", sep="\t", index=False)

def run_one(name, ad):
    try:
        li.mt.rank_aggregate(ad, groupby="grp", resource=sub_res, expr_prop=0.1, use_raw=False, n_perms=100, verbose=False, seed=0)
        x = ad.uns["liana_res"][["source", "target", "ligand_complex", "receptor_complex", "magnitude_rank"]].copy()
        for c in ["source", "target", "ligand_complex", "receptor_complex"]: x[c] = x[c].astype(str)
        x["sample"] = name; x["groups"] = ",".join(sorted(ad.obs["grp"].unique()))
        return x
    except Exception as e:
        return pd.DataFrame({"sample": [name], "err": [repr(e)[:200]]})

def prep(a, samples, labcol, tag):
    rng = np.random.default_rng(0); out = []
    for s in samples:
        m = (a.obs["sample"].astype(str) == s).values
        sub = a[m]
        grp = np.array([""] * sub.n_obs, dtype=object)
        lab = sub.obs[labcol].astype(str).values
        grp[np.isin(lab, EPI)] = "Epithelial"
        for L, labs in LIN.items(): grp[np.isin(lab, labs)] = L
        keep = []
        for g in np.unique(grp[grp != ""]):
            idx = np.where(grp == g)[0]
            if len(idx) >= 20: keep.append(rng.choice(idx, min(len(idx), 500), replace=False))
        if not keep: continue
        keep = np.sort(np.concatenate(keep)); ad = sub[keep].copy(); ad.obs["grp"] = pd.Categorical(grp[keep]); ad.obs = ad.obs[["grp"]]
        ad.layers = {}; ad.obsm = {}; ad.varm = {}; ad.obsp = {}
        out.append((s, ad))
    return out

def score_table(runs, samples):
    df = pd.concat([r for r in runs if "magnitude_rank" in r.columns], ignore_index=True)
    df["score"] = 1 - df["magnitude_rank"]
    key = pairs.rename(columns={"source": "src", "target": "tgt"})
    df = df.merge(key.rename(columns={"ligand": "ligand_complex", "receptor": "receptor_complex"}), left_on=["ligand_complex", "receptor_complex", "source", "target"],
                  right_on=["ligand_complex", "receptor_complex", "src", "tgt"]).drop(columns=["src", "tgt"])
    df["pair"] = df["ligand_complex"] + "→" + df["receptor_complex"] + "|" + df["source"] + ">" + df["target"]
    return df

def matrix(df, runs):
    groups = {r["sample"].iloc[0]: set(r["groups"].iloc[0].split(",")) for r in runs if "groups" in r.columns}
    M = df.pivot_table(index="pair", columns="sample", values="score", aggfunc="max")
    allp = pairs.assign(pair=lambda x: x.ligand + "→" + x.receptor + "|" + x.source + ">" + x.target)
    M = M.reindex(allp["pair"].unique())
    for s, gs in groups.items():
        for p in M.index:
            src, tgt = p.split("|")[1].split(">")
            if src in gs and tgt in gs and pd.isna(M.loc[p, s]): M.loc[p, s] = 0.0
    return M

def main():
    log("===== 开始 =====")
    a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
    a.obs["sample"] = a.obs["sample"].astype(str)
    ep = np.asarray(a[:, "EPCAM"].layers["counts"].todense()).ravel()
    lab = a.obs["final_label"].astype(str)
    a = a[~(~lab.isin(EPI) & (ep >= 2)).values].copy()
    meta = a.obs.groupby("sample").agg(patient=("patient", "first"), stage=("stage", "first"))
    meta = meta.astype(str)
    meta = meta[meta.stage != "Normal"]
    log(f"Peng 样本 {len(meta)}")
    jobs = prep(a, list(meta.index), "final_label", "peng"); del a
    runs = Parallel(n_jobs=5)(delayed(run_one)(n, ad) for n, ad in jobs)
    err = [r for r in runs if "err" in r.columns]
    log(f"Peng 运行完成 {len(runs)-len(err)}/{len(runs)}；错误 {len(err)} {err[0]['err'].iloc[0] if err else ''}")
    df = score_table(runs, meta.index); M = matrix(df, runs)
    M.to_csv(f"{OUT}/peng_pair_scores.tsv", sep="\t")
    pat = meta.loc[M.columns, "patient"]; grp = np.where(meta.loc[M.columns, "stage"] == "LUAD", "LUAD", "Precursor")
    rows = []
    for p in M.index:
        d = {}
        for pt in pat.unique():
            cols = M.columns[(pat == pt).values]; g = grp[(pat == pt).values]
            a1 = M.loc[p, cols[g == "Precursor"]].mean(); a2 = M.loc[p, cols[g == "LUAD"]].mean()
            if pd.notna(a1) and pd.notna(a2): d[pt] = a2 - a1
        d = pd.Series(d)
        if len(d) >= 8 and (d != 0).sum() >= 5:
            rows.append(dict(pair=p, n_pairs=len(d), mean_diff=d.mean(), p=wilcoxon(d).pvalue, n_up=int((d > 0).sum())))
    t = pd.DataFrame(rows)
    if len(t):
        t["q"] = multipletests(t["p"], method="fdr_bh")[1]
        t["sig_up"] = (t["q"] < 0.05) & (t["mean_diff"] > 0)
    t.to_csv(f"{OUT}/peng_paired_test.tsv", sep="\t", index=False)
    log(f"Peng：可检验相互作用 {len(t)}；q<0.05 且 LUAD 更高 {int(t['sig_up'].sum()) if len(t) else 0}；q<0.05 且前体更高 {int(((t['q']<0.05)&(t['mean_diff']<0)).sum()) if len(t) else 0}")

    # ---------- GSE164789 验证 ----------
    v = sc.read_h5ad(f"{P}/results/05_validation/GSE164789/GSE164789_annotated.h5ad")
    ann = pd.read_csv(f"{P}/results/05_validation/GSE164789/cluster_annotation_auto.tsv", sep="\t", index_col=0)
    v.obs["fl"] = v.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
    v.obs["sample"] = v.obs["sample"].astype(str)
    ep = np.asarray(v[:, "EPCAM"].layers["counts"].todense()).ravel()
    v = v[(v.obs["tissue"] == "T").values & ~(~v.obs["fl"].isin(EPI).values & (ep >= 2))].copy()
    vm = v.obs.groupby("sample").agg(patient=("patient", "first"), histology=("histology", "first"))
    vm = vm.astype(str)
    jobs = prep(v, list(vm.index), "fl", "gse164789"); del v
    runs = Parallel(n_jobs=5)(delayed(run_one)(n, ad) for n, ad in jobs)
    err = [r for r in runs if "err" in r.columns]; log(f"GSE164789 运行完成 {len(runs)-len(err)}/{len(runs)}；错误 {len(err)}")
    df = score_table(runs, vm.index); M2 = matrix(df, runs); M2.to_csv(f"{OUT}/gse164789_pair_scores.tsv", sep="\t")
    unit = (vm.loc[M2.columns, "patient"] + "_" + vm.loc[M2.columns, "histology"]).values
    M3 = M2.T.groupby(unit).mean().T
    cond = pd.Series(["IAC" if u.endswith("_IAC") else "Precursor" for u in M3.columns], index=M3.columns)
    sig = t[t["sig_up"]]["pair"].tolist() if len(t) else []
    rows = []
    for p in sig:
        x = M3.loc[p]; i1 = x[cond == "IAC"].dropna(); i0 = x[cond == "Precursor"].dropna()
        if len(i1) >= 3 and len(i0) >= 3:
            rows.append(dict(pair=p, n_IAC=len(i1), n_pre=len(i0), mean_diff=i1.mean() - i0.mean(), p1=mannwhitneyu(i1, i0, alternative="greater").pvalue))
    vt = pd.DataFrame(rows)
    if len(vt):
        vt["q"] = multipletests(vt["p1"], method="fdr_bh")[1]; vt["pass"] = (vt["q"] < 0.05) & (vt["mean_diff"] > 0)
        vt = vt.merge(t[["pair", "mean_diff", "q"]].rename(columns={"mean_diff": "peng_diff", "q": "peng_q"}), on="pair")
    vt.to_csv(f"{OUT}/gse164789_validation.tsv", sep="\t", index=False)
    log(f"GSE164789 验证：Peng 显著上升 {len(sig)}；可检验 {len(vt)}；方向一致 {int((vt['mean_diff']>0).sum()) if len(vt) else 0}；q<0.05 通过 {int(vt['pass'].sum()) if len(vt) else 0}")
    if len(vt):
        log("通过者：" + ", ".join(vt[vt['pass']].sort_values('q')['pair'].head(40)))
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
