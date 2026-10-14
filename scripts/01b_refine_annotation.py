#!/usr/bin/env python
"""
第一步（下）：在已确认的 final_label（cluster_annotation_auto.tsv）上继续精修注释。

用法：
    cd ~/全部数据库/GGN_target && nohup /opt/anaconda3/envs/ggn_target/bin/python scripts/01b_refine_annotation.py > results/01b_console.txt 2>&1 &

处理顺序（剔除原因按此顺序归属，一个细胞只记第一个原因）：
  A. 剔除已确认的双细胞簇（3）和血小板簇（17）：final_label = Doublet/ambiguous、Platelet/Megakaryocyte
  B. doublet_score > 0.25 的细胞（选择依据见下）
  C. 拆分原簇 5（增殖混合簇）和原簇 19（基质-上皮混合簇）：子集内用 X_pca_harmony 重建邻居图，
     leiden 分辨率 0.5，每个亚簇用 01 脚本的 MARKERS 打分，按固定规则定类；混合或无明确标志的亚簇剔除
  D. 跨谱系共表达：上皮细胞中 PTPRC 原始计数 >= 2；免疫细胞中 EPCAM 原始计数 >= 2

doublet_score 阈值依据（运行前先查看了分布）：全体 99 百分位 0.227、99.9 百分位 0.385；
0.25 约为 99.3 百分位，剔除全体的 0.74%，且各样本 99 百分位多在 0.24–0.34，分布为单峰带右尾，
没有清晰的双峰，故采用接近全局 99 百分位的 0.25（保守，宁可少剔）。图见 doublet_score_by_sample.png。

亚簇定类规则（预先写死；2026-10-08 第一次试跑发现仅靠标志基因得分会被环境 RNA 抬高上皮得分，故加入 EPCAM/PTPRC 均值）：
  对每个亚簇计算：各谱系组内最高标志基因得分；EPCAM 均值；PTPRC 均值。
  1) 上皮证据 = EPCAM 均值 >= 1.0，或（EPCAM 均值 >= 0.8 且 Epithelial 得分 >= 0.3）；免疫证据 = PTPRC 均值 >= 0.4。
  2) 上皮证据与免疫证据同时存在 → mixed，剔除。
  3) 仅有上皮证据 → 保留：原簇 5 标 Proliferating epithelial；原簇 19 取上皮内得分最高的类型。
  4) 仅有免疫证据 → 在非上皮谱系中取得分最高者（需 >= 0.3，否则 no_clear_marker 剔除）；
     第二名 >= 第一名的 50% → mixed 剔除；原簇 5 加 Proliferating 前缀。
  5) 两种证据都没有 → 在非上皮谱系中判断：无 >= 0.3 者 → no_clear_marker 剔除；第二名 >= 第一名 50% → mixed 剔除；
     EPCAM 均值 >= 0.4 且 Epithelial 得分 >= 第一名 50%（上皮污染的基质/内皮）→ mixed 剔除；否则归入最高谱系。
  Plasma 的标志基因去掉 XBP1（上皮细胞也表达）。
跨谱系共表达（步骤 D）：上皮细胞 PTPRC 原始计数 >= 2；免疫细胞 EPCAM 原始计数 >= 2。
  SFTPB 在 43.6% 的免疫细胞中计数 >= 2（环境 RNA），不适合作为跨谱系标志，已弃用，仅在日志中记录。
"""
import os, json, time, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

BASE = os.path.expanduser("~/全部数据库/GGN_target")
OUT = os.path.join(BASE, "results/01_snRNA")
SEED = 0
CL = "leiden_1.0"
DOUBLET_THR = 0.25
STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]
SPLIT = {"5": "Proliferating", "19": "Mixed-stromal-epithelial"}   # 原簇号 → 说明

LOG = open(os.path.join(OUT, "01b_log.txt"), "a", encoding="utf-8")
T0 = time.time()
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {msg}"
    print(line, flush=True); LOG.write(line + "\n"); LOG.flush()

MARKERS = {
    "AT2":            ["SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"],
    "AT1":            ["AGER", "RTKN2", "CLIC5", "EMP2"],
    "Club/secretory": ["SCGB1A1", "SCGB3A2", "SCGB3A1", "BPIFB1"],
    "Ciliated":       ["FOXJ1", "CAPS", "TPPP3", "RSPH1"],
    "Basal":          ["KRT5", "TP63", "KRT15"],
    "Endothelial":    ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Lymphatic EC":   ["PROX1", "CCL21", "MMRN1"],
    "Fibroblast":     ["COL1A1", "COL1A2", "DCN", "LUM"],
    "SMC/pericyte":   ["ACTA2", "MYH11", "RGS5", "PDGFRB"],
    "Mesothelial":    ["MSLN", "UPK3B", "WT1"],
    "T cell":         ["CD3E", "CD3D", "TRAC", "CD2"],
    "NK":             ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "B cell":         ["MS4A1", "CD79A", "BANK1"],
    "Plasma":         ["JCHAIN", "MZB1", "IGHG1"],
    "Macrophage":     ["CD68", "C1QA", "C1QB", "MARCO", "APOE"],
    "Monocyte":       ["CD14", "VCAN", "S100A8"],
    "cDC":            ["CLEC9A", "CD1C", "FCER1A", "CLEC10A"],
    "pDC":            ["LILRA4", "IL3RA", "CLEC4C"],
    "Mast":           ["CPA3", "KIT", "MS4A2", "HDC"],
    "Neutrophil":     ["CSF3R", "FCGR3B", "CXCR2"],
}
LINEAGE = {
    "Epithelial": ["AT2", "AT1", "Club/secretory", "Ciliated", "Basal"],
    "T/NK": ["T cell", "NK"],
    "Myeloid": ["Macrophage", "Monocyte", "cDC", "Neutrophil", "Mast"],
    "B/Plasma": ["B cell", "Plasma", "pDC"],
    "Stromal": ["Fibroblast", "SMC/pericyte", "Mesothelial"],
    "Endothelial": ["Endothelial", "Lymphatic EC"],
}
PROLIF_NAME = {"Epithelial": "Proliferating epithelial", "T/NK": "Proliferating T/NK",
               "Myeloid": "Proliferating myeloid", "B/Plasma": "Proliferating B/Plasma",
               "Stromal": "Proliferating stromal", "Endothelial": "Proliferating endothelial"}
EPI_LABELS = ["AT2", "AT1", "Ciliated", "Basal", "Proliferating epithelial",
              "AT2/Club-like epithelial", "AT2-like (patient-specific)"]
STROMAL = ["Fibroblast", "SMC/pericyte", "Mesothelial"]
ENDO = ["Endothelial", "Lymphatic EC"]
DROP_LABELS = {"Doublet/ambiguous": "cluster_doublet", "Platelet/Megakaryocyte": "cluster_platelet"}
PRESENT_THR, MIXED_RATIO = 0.3, 0.5

def compartment(l):
    if l in EPI_LABELS: return "Epithelial"
    if l in STROMAL or l == "Proliferating stromal": return "Stromal"
    if l in ENDO or l == "Proliferating endothelial": return "Endothelial"
    return "Immune"

def counts_of(adata, genes, mask=None):
    g = [x for x in genes if x in adata.var_names]
    X = adata[:, g].layers["counts"]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return pd.DataFrame(X, index=adata.obs_names, columns=g)

def main():
    log("===== 开始 =====")
    log(f"scanpy {sc.__version__}")
    adata = sc.read_h5ad(os.path.join(OUT, "snRNA_annotated.h5ad"))
    log(f"读取：{adata.shape}；obs 列：{list(adata.obs.columns)}")
    ann = pd.read_csv(os.path.join(OUT, "cluster_annotation_auto.tsv"), sep="\t", index_col=0)
    ann.index = ann.index.astype(str)
    assert (ann["final_label"].fillna("") != "").all()
    obs = adata.obs
    obs["orig_cluster"] = obs[CL].astype(str)
    obs["final_label"] = obs["orig_cluster"].map(ann["final_label"]).astype(str)
    obs["remove_reason"] = ""

    # ---------- A. 已确认的剔除簇 ----------
    for l, r in DROP_LABELS.items():
        obs.loc[(obs["final_label"] == l) & (obs["remove_reason"] == ""), "remove_reason"] = r

    # ---------- B. doublet_score ----------
    samples = [s for s in obs["sample"].cat.categories] if hasattr(obs["sample"], "cat") else sorted(obs["sample"].unique())
    fig, ax = plt.subplots(figsize=(18, 4.5))
    ax.boxplot([obs.loc[obs["sample"] == s, "doublet_score"].values for s in samples], showfliers=False)
    ax.set_xticks(range(1, len(samples) + 1)); ax.set_xticklabels(samples, rotation=90, fontsize=6)
    ax.axhline(DOUBLET_THR, color="r", ls="--"); ax.set_ylabel("doublet_score")
    fig.savefig(os.path.join(OUT, "doublet_score_by_sample.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
    m = (obs["doublet_score"] > DOUBLET_THR) & (obs["remove_reason"] == "")
    obs.loc[m, "remove_reason"] = "doublet_score"
    log(f"doublet_score > {DOUBLET_THR}：{int(m.sum())} 个细胞（全体 {m.mean()*100:.2f}%）")

    # ---------- C. 拆分原簇 5 与 19 ----------
    sub_rows, sub_marker_rows = [], []
    for oc, why in SPLIT.items():
        keep = (obs["orig_cluster"] == oc) & (obs["remove_reason"] == "")
        sub = adata[keep.values].copy()
        log(f"原簇 {oc}（{why}）：{sub.n_obs} 个细胞，重建邻居图并 leiden 0.5")
        sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=30, random_state=SEED)
        sc.tl.leiden(sub, resolution=0.5, key_added="sub", flavor="igraph", n_iterations=2,
                     directed=False, random_state=SEED)
        for k, v in MARKERS.items():
            gs = [g for g in v if g in sub.var_names]
            sc.tl.score_genes(sub, gs, score_name=f"score::{k}", random_state=SEED)
        sc_cols = [f"score::{k}" for k in MARKERS]
        cs = sub.obs.groupby("sub", observed=True)[sc_cols].mean()
        cs.columns = list(MARKERS)
        gen = [g for g in ["EPCAM", "KRT8", "KRT19", "PTPRC", "MKI67", "TOP2A"] if g in sub.var_names]
        gx = sub[:, gen].X
        gx = gx.toarray() if sp.issparse(gx) else np.asarray(gx)
        gm = pd.DataFrame(gx, index=sub.obs_names, columns=gen).groupby(sub.obs["sub"].values).mean()
        sc.tl.rank_genes_groups(sub, "sub", method="wilcoxon", n_genes=10)
        topg = {c: ",".join(sub.uns["rank_genes_groups"]["names"][c][:10])
                for c in sub.uns["rank_genes_groups"]["names"].dtype.names}
        decision = {}
        for c in cs.index:
            lin = {L: cs.loc[c, ts].max() for L, ts in LINEAGE.items()}
            epcam, ptprc = float(gm.loc[c, "EPCAM"]), float(gm.loc[c, "PTPRC"])
            epi_ev = epcam >= 1.0 or (epcam >= 0.8 and lin["Epithelial"] >= PRESENT_THR)
            imm_ev = ptprc >= 0.4
            non = sorted([L for L in lin if L != "Epithelial"], key=lin.get, reverse=True)
            top, sec = non[0], non[1]
            second_close = lin[sec] >= MIXED_RATIO * lin[top]
            if epi_ev and imm_ev:
                dec, lab = "remove", "mixed"
            elif epi_ev:
                t = max(LINEAGE["Epithelial"], key=lambda x: cs.loc[c, x])
                dec, lab = "keep", ("Proliferating epithelial" if oc == "5" else t)
            elif imm_ev:
                if lin[top] < PRESENT_THR:
                    dec, lab = "remove", "no_clear_marker"
                elif second_close:
                    dec, lab = "remove", "mixed"
                else:
                    t = max(LINEAGE[top], key=lambda x: cs.loc[c, x])
                    dec, lab = "keep", (PROLIF_NAME[top] if oc == "5" else t)
            else:
                if lin[top] < PRESENT_THR:
                    dec, lab = "remove", "no_clear_marker"
                elif second_close or (epcam >= 0.4 and lin["Epithelial"] >= MIXED_RATIO * lin[top]):
                    dec, lab = "remove", "mixed"
                else:
                    t = max(LINEAGE[top], key=lambda x: cs.loc[c, x])
                    dec, lab = "keep", (PROLIF_NAME[top] if oc == "5" else t)
            order = sorted(lin, key=lin.get, reverse=True)
            decision[c] = (dec, lab)
            sub_rows.append(dict(orig_cluster=oc, sub=c, n_cells=int((sub.obs["sub"] == c).sum()),
                                 n_patients=sub.obs.loc[sub.obs["sub"] == c, "patient"].nunique(),
                                 decision=dec, label=lab,
                                 lineage_scores="; ".join(f"{L}={lin[L]:.2f}" for L in order),
                                 **{f"mean_{g}": round(float(gm.loc[c, g]), 3) for g in gen},
                                 top_markers=topg[c]))
        for c, (dec, lab) in decision.items():
            idx = sub.obs_names[sub.obs["sub"] == c]
            if dec == "keep":
                obs.loc[idx, "final_label"] = lab
            else:
                obs.loc[idx, "final_label"] = f"REMOVED_{lab}"
                obs.loc[idx, "remove_reason"] = f"subcluster_{lab}"
            log(f"  原簇 {oc} 亚簇 {c}：{len(idx)} 个细胞 → {dec} / {lab}")
        del sub
    pd.DataFrame(sub_rows).to_csv(os.path.join(OUT, "subcluster_5_19_decisions.tsv"), sep="\t", index=False)

    # ---------- D. 跨谱系共表达 ----------
    alive = obs["remove_reason"] == ""
    is_epi = obs["final_label"].isin(EPI_LABELS) & alive
    is_imm = obs["final_label"].map(compartment).eq("Immune") & alive
    cp = counts_of(adata, ["PTPRC", "EPCAM", "SFTPB"])
    m1 = is_epi & (cp["PTPRC"] >= 2)
    m2 = is_imm & (cp["EPCAM"] >= 2)
    log(f"（仅记录）免疫细胞 SFTPB>=2 的比例 {(is_imm & (cp['SFTPB'] >= 2)).sum()/is_imm.sum()*100:.1f}%，为环境 RNA，不用于剔除")
    for nm, mm, den in [("上皮 PTPRC", m1, is_epi), ("免疫 EPCAM", m2, is_imm)]:   # 合理性上限：上皮-免疫双细胞不应超过 5%
        if mm.sum() / den.sum() > 0.05:
            raise RuntimeError(f"跨谱系规则「{nm}」剔除比例 {mm.sum()/den.sum()*100:.1f}% > 5%，规则可能有问题，停止检查")
    obs.loc[m1, "remove_reason"] = "cross_epi_PTPRC"
    obs.loc[m2, "remove_reason"] = "cross_imm_EPCAM"
    log(f"跨谱系：上皮含 PTPRC>=2 {int(m1.sum())}（占上皮 {m1.sum()/is_epi.sum()*100:.1f}%）；"
        f"免疫含 EPCAM>=2 {int(m2.sum())}（占免疫 {m2.sum()/is_imm.sum()*100:.1f}%）")

    # ---------- 剔除汇总（先报告，再执行） ----------
    reasons = ["cluster_doublet", "cluster_platelet", "doublet_score", "subcluster_mixed",
               "subcluster_no_clear_marker", "cross_epi_PTPRC", "cross_imm_EPCAM"]
    tab = pd.crosstab(obs["sample"], obs["remove_reason"]).drop(columns="", errors="ignore")
    for r in reasons:
        if r not in tab.columns: tab[r] = 0
    tab = tab[reasons]
    tab["n_before"] = obs.groupby("sample", observed=True).size()
    tab["n_removed"] = tab[reasons].sum(axis=1)
    tab["frac_removed"] = (tab["n_removed"] / tab["n_before"]).round(4)
    tab.to_csv(os.path.join(OUT, "removed_cells_summary.tsv"), sep="\t")
    log("各原因合计剔除：" + json.dumps({r: int(tab[r].sum()) for r in reasons}, ensure_ascii=False))
    log(f"合计剔除 {int(tab['n_removed'].sum())} / {len(obs)}（{tab['n_removed'].sum()/len(obs)*100:.2f}%）")
    high = tab[tab["frac_removed"] > 0.15]
    log("剔除比例 >15% 的样本：" + (", ".join(f"{i}({r*100:.1f}%)" for i, r in high["frac_removed"].items()) or "无"))
    for s, r in tab.iterrows():
        log(f"  {s}: 剔除 {int(r['n_removed'])}/{int(r['n_before'])}（{r['frac_removed']*100:.1f}%）")

    # ---------- 执行剔除并写出 ----------
    adata = adata[(obs["remove_reason"] == "").values].copy()
    adata.obs["final_label"] = pd.Categorical(adata.obs["final_label"].astype(str))
    adata.obs["final_compartment"] = pd.Categorical(adata.obs["final_label"].astype(str).map(compartment))
    log(f"剔除后：{adata.n_obs} 个细胞；final_label：{adata.obs['final_label'].value_counts().to_dict()}")

    lab = adata.obs.groupby("final_label", observed=True).agg(
        n_cells=("sample", "size"), n_samples=("sample", "nunique"), n_patients=("patient", "nunique"),
        origin_clusters=("orig_cluster", lambda x: ",".join(sorted(x.unique(), key=int))),
        compartment=("final_compartment", "first"))
    lab.to_csv(os.path.join(OUT, "cluster_annotation_final.tsv"), sep="\t")

    comp = pd.crosstab(adata.obs["sample"], adata.obs["final_label"])
    comp.to_csv(os.path.join(OUT, "composition_final_counts.tsv"), sep="\t")
    comp.div(comp.sum(axis=1), axis=0).round(4).to_csv(os.path.join(OUT, "composition_final_frac.tsv"), sep="\t")

    fig = sc.pl.umap(adata, color="final_label", show=False, return_fig=True, frameon=False)
    fig.savefig(os.path.join(OUT, "umap_final_label.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
    genes_dot = {k: [g for g in v if g in adata.var_names][:2] for k, v in MARKERS.items()}
    dp = sc.pl.dotplot(adata, genes_dot, groupby="final_label", standard_scale="var", show=False, return_fig=True)
    dp.savefig(os.path.join(OUT, "dotplot_final.png"), dpi=130, bbox_inches="tight"); plt.close("all")

    adata.obs = adata.obs.drop(columns=["remove_reason"])
    adata.write_h5ad(os.path.join(OUT, "snRNA_final.h5ad"), compression="lzf")
    epi = adata[(adata.obs["final_compartment"] == "Epithelial").values].copy()
    epi.write_h5ad(os.path.join(OUT, "epithelial_cells.h5ad"), compression="lzf")
    ps = epi.obs.groupby(["sample"], observed=True).agg(
        patient=("patient", "first"), stage=("stage", "first"), n_epi=("sample", "size"))
    ps.to_csv(os.path.join(OUT, "epithelial_cells_per_sample.tsv"), sep="\t")
    log(f"上皮细胞 {epi.n_obs} 个；每样本上皮细胞数 <50 的样本：{list(ps.index[ps['n_epi'] < 50]) or '无'}")
    log("===== 完成 =====")

if __name__ == "__main__":
    main()
