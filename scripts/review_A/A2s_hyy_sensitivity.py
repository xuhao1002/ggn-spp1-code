"""A2s（偏离记录 D2，post hoc）GSE164789 一名患者（HYY）贡献两个单位的敏感性分析。环境 ggn_target。
管线与 14_microenv_layer1.py 的外部验证部分逐项一致，只改变纳入的单位：
  V0  全部 30 单位（复现原结果：SPP1 log2FC 2.20，单侧 P 0.0014）
  V1  去掉 HYY 两个单位（14 IAC 对 14 前体）
  V2a 只保留 HYY 的 IAC 单位（15 对 14）
  V2b 只保留 HYY 的 MIA 单位（14 对 15）
判定（运行前写定于 D2）：V1、V2a、V2b 三者全部"同向且 BH q < 0.05" → 稳健；任一不满足 → 不稳健。
主结果（V0）在两种情况下都不替换。"""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import numpy as np, pandas as pd, scanpy as sc
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A2")
DS = "GSE164789"; HYY = "HYY"; MIN_CELLS = 50
MYE = ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"]
ORIG = {"lfc": 2.202066220713526, "p1": 0.001443149362762689, "q": 0.01754675355012904}

# ---- 与 14_microenv_layer1.py 完全相同的两个函数 ----
def filt(c):
    n = len(c); k = ((c >= 10).sum(axis=0) >= np.ceil(n / 3)) & (c.sum(axis=0) >= 50); return c.loc[:, k.values]

def deseq(counts, meta, design, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    c = filt(counts)
    dds = DeseqDataSet(counts=c, metadata=meta, design=design, refit_cooks=True, n_cpus=6, quiet=True); dds.deseq2()
    ds = DeseqStats(dds, contrast=contrast, quiet=True); ds.summary()
    r = ds.results_df.copy(); r.index.name = "gene"; return r

# ---- BH 校正范围：14 中髓系、发现队列显著且通过来源过滤的候选基因（与原分析相同） ----
cand = pd.read_csv(f"{P}/results/14_microenv/microenv_layer1_candidates.tsv", sep="\t", index_col=0)
cand = cand[cand["lineage"] == "Myeloid"]
assert "SPP1" in cand.index
disc_sign = np.sign(cand["log2FoldChange"])

# ---- 构建单位（同 14） ----
d = f"{P}/results/05_validation/{DS}"
v = sc.read_h5ad(f"{d}/{DS}_annotated.h5ad"); log("载入", v.shape)
ann = pd.read_csv(f"{d}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
fl = v.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str).values
ep = np.asarray(v[:, "EPCAM"].layers["counts"].todense()).ravel()
keep = np.isin(fl, MYE) & (ep < 2) & (v.obs["tissue"].astype(str).values == "T")
v = v[keep].copy()
v.obs["unit"] = v.obs["patient"].astype(str) + "_" + v.obs["histology"].astype(str)
v.obs["condition"] = np.where(v.obs["histology"].astype(str) == "IAC", "IAC", "Precursor")
X = v.layers["counts"].tocsr()
groups, meta = {}, []
for u, idx in v.obs.groupby("unit", observed=True).indices.items():
    if len(idx) < MIN_CELLS: continue
    groups[u] = idx
    meta.append(dict(unit=u, patient=str(v.obs["patient"].iloc[idx[0]]),
                     condition=v.obs["condition"].iloc[idx[0]], n_cells=len(idx)))
meta = pd.DataFrame(meta).set_index("unit")
counts_all = pd.DataFrame({k: np.asarray(X[i].sum(axis=0)).ravel() for k, i in groups.items()},
                          index=v.var_names).T.astype(int)
del v, X
hyy = meta.index[meta.patient == HYY].tolist(); log("HYY 单位", hyy)
hyy_iac = [u for u in hyy if meta.loc[u, "condition"] == "IAC"]
hyy_pre = [u for u in hyy if meta.loc[u, "condition"] == "Precursor"]
assert len(hyy_iac) == 1 and len(hyy_pre) == 1, hyy

VARIANTS = {"V0_all30": [], "V1_drop_both_HYY": hyy, "V2a_keep_HYY_IAC": hyy_pre, "V2b_keep_HYY_MIA": hyy_iac}
rows = []
for name, drop in VARIANTS.items():
    m = meta.drop(index=drop).copy()
    c = counts_all.loc[m.index]
    m["condition"] = pd.Categorical(m["condition"], categories=["Precursor", "IAC"])
    r = deseq(c, m[["condition"]], "~ condition", ["condition", "IAC", "Precursor"])
    vv = r.reindex(cand.index)
    same = np.sign(vv["log2FoldChange"]) == disc_sign
    p1 = np.where(same, vv["pvalue"] / 2, 1 - vv["pvalue"] / 2)
    ok = vv["pvalue"].notna().values
    q = np.full(len(cand), np.nan)
    if ok.sum() > 0: q[ok] = multipletests(p1[ok], method="fdr_bh")[1]
    i = list(cand.index).index("SPP1")
    row = dict(variant=name, n_IAC=int((m.condition == "IAC").sum()), n_Precursor=int((m.condition == "Precursor").sum()),
               n_patients=int(m.patient.nunique()), n_units=int(len(m)),
               SPP1_log2FC=float(vv["log2FoldChange"].iloc[i]), SPP1_lfcSE=float(vv["lfcSE"].iloc[i]),
               SPP1_p_two_sided=float(vv["pvalue"].iloc[i]), SPP1_p_one_sided=float(p1[i]), SPP1_q_BH=float(q[i]),
               SPP1_same_direction=bool(same.iloc[i]), n_genes_in_BH=int(ok.sum()))
    row["passes_replication_criterion"] = bool(row["SPP1_same_direction"] and row["SPP1_q_BH"] < 0.05)
    rows.append(row); log(name, {k: row[k] for k in ["n_IAC", "n_Precursor", "SPP1_log2FC", "SPP1_p_one_sided", "SPP1_q_BH"]})

R = pd.DataFrame(rows); R.to_csv(OUT + "A2s_HYY_sensitivity.tsv", sep="\t", index=False)
v0 = R.set_index("variant").loc["V0_all30"]
repro = bool(abs(v0.SPP1_log2FC - ORIG["lfc"]) < 0.01 and abs(v0.SPP1_p_one_sided - ORIG["p1"]) / ORIG["p1"] < 0.05)
sens = R[R.variant != "V0_all30"]
robust = bool(sens.passes_replication_criterion.all())
dump({"plan_sha256": sha, "deviation": "D2 (post hoc)",
      "original_reported": ORIG, "V0_reproduces_original": repro,
      "results": rows,
      "decision_rule": "V1、V2a、V2b 全部同向且 BH q<0.05 → 稳健；任一不满足 → 不稳健（运行前写定）",
      "robust": robust if repro else None,
      "note": ("V0 未能复现原结果：先查原因，不解读 V1–V2。" if not repro else
               ("三个敏感性版本全部满足复制标准：采用'稳健'版措辞。" if robust else
                "至少一个敏感性版本不满足复制标准：采用'不稳健'版措辞，主 P 值不换，结果与摘要同步退。"))},
     OUT + "A2s_HYY_sensitivity.json")
log("A2s 完成")
