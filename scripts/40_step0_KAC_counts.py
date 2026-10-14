#!/usr/bin/env python
"""第四部分第 0 步（只读核查，不做任何统计检验）：每名患者前体侧与 Normal 侧的 KAC 细胞数、AT2L-other 细胞数，以及两侧都 >=50 的可配对数。
- KAC 细胞 = AT2 谱系（AT2、AT2/Club-like、增殖上皮）细胞中 KAC 得分 > 阈值；阈值直接读取第二阶段已固定的数值（results/20_stage2/discovery/run_log.txt），不重新计算。
- KAC 得分用与第二阶段完全相同的做法重现（scanpy score_genes，Table S3 人类 KAC 91/102 个可用基因，random_state=0，作用于全部 AT2 谱系细胞）；
  仅用于计数；并用第二阶段的 KAC 比例表（cell_fractions_AT2L.tsv）核对重现是否一致。
- 前体 = 同一患者 AAH/AIS/MIA 合并；Normal = Normal + Normal1。另给出 AAH 与 AIS/MIA 分开的计数（描述性）。
输出 results/40_drivers/state_matched/step0_KAC_counts.tsv、step0_summary.txt；日志 results/40_drivers/run_log.txt"""
import os, re, time, warnings
import numpy as np, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")
P = os.path.expanduser("~/全部数据库/GGN_target"); OUT = f"{P}/results/40_drivers"
LOG = open(f"{OUT}/run_log.txt", "a", encoding="utf-8"); T0 = time.time()
def log(m):
    line = f"[{time.strftime('%H:%M:%S')} +{(time.time()-T0)/60:5.1f} min] {m}"; print(line, flush=True); LOG.write(line + "\n"); LOG.flush()
txt = open(f"{P}/results/20_stage2/discovery/run_log.txt", encoding="utf-8").read()
m = re.search(r"阈值（Normal AT2 谱系细胞第 90 百分位）：\{'KAC': np.float64\(([0-9.]+)\), 'CAS': np.float64\(([0-9.]+)\)\}", txt)
THR_KAC, THR_CAS = float(m.group(1)), float(m.group(2))
log(f"===== 第四部分 第 0 步（只读核查）；读取第二阶段固定阈值：KAC={THR_KAC}，CAS={THR_CAS} =====")
AT2L = ["AT2", "AT2/Club-like epithelial", "Proliferating epithelial"]
kac = [l.strip() for l in open(f"{P}/results/20_stage2/peng_TableS3_KAC_human.txt")]
a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
o = a.obs
m_at = o["final_label"].astype(str).isin(AT2L).values
s2 = a[m_at].copy(); s2.obs = s2.obs[["patient", "stage"]].copy(); s2.obs["patient"] = s2.obs["patient"].astype(str); s2.obs["stage"] = s2.obs["stage"].astype(str)
del a
pres = [g for g in kac if g in s2.var_names]
sc.tl.score_genes(s2, pres, score_name="KAC", random_state=0)
sc_ = s2.obs["KAC"].values
# 重现性核对（不重算阈值）：Normal 细胞得分第 90 百分位应与固定阈值一致
rec = np.percentile(sc_[(s2.obs["stage"] == "Normal").values], 90)
log(f"重现性核对：KAC 基因 {len(pres)} 个；Normal AT2 谱系细胞得分的第 90 百分位（仅核对）= {rec:.10f}，与固定阈值 {THR_KAC:.10f} 的差 {abs(rec-THR_KAC):.2e}")
s2.obs["is_kac"] = sc_ > THR_KAC
grp = s2.obs["stage"].map(lambda s: "Normal" if s == "Normal" else ("Precursor" if s in ("AAH", "AIS", "MIA") else "LUAD"))
s2.obs["grp"] = grp.values
rows = []
for p in sorted(s2.obs["patient"].unique(), key=lambda x: int(x[1:])):
    sub = s2.obs[s2.obs["patient"] == p]
    def cnt(mask): return int(mask.sum())
    nrm = sub["grp"] == "Normal"; pre = sub["grp"] == "Precursor"
    aah = sub["stage"] == "AAH"; am = sub["stage"].isin(["AIS", "MIA"])
    rows.append(dict(patient=p, normal_AT2L=cnt(nrm), normal_KAC=cnt(nrm & sub.is_kac), normal_other=cnt(nrm & ~sub.is_kac),
                     precursor_AT2L=cnt(pre), precursor_KAC=cnt(pre & sub.is_kac), precursor_other=cnt(pre & ~sub.is_kac),
                     AAH_KAC=cnt(aah & sub.is_kac), AAH_other=cnt(aah & ~sub.is_kac), AISMIA_KAC=cnt(am & sub.is_kac), AISMIA_other=cnt(am & ~sub.is_kac)))
t = pd.DataFrame(rows).set_index("patient")
t["pair_KAC"] = (t.normal_KAC >= 50) & (t.precursor_KAC >= 50)
t["pair_other"] = (t.normal_other >= 50) & (t.precursor_other >= 50)
t.to_csv(f"{OUT}/state_matched/step0_KAC_counts.tsv", sep="\t")
# 与第二阶段 KAC 比例表核对
fr = pd.read_csv(f"{P}/results/20_stage2/discovery/cell_fractions_AT2L.tsv", sep="\t")
chk = []
for _, r in fr.iterrows():
    p, g = r["patient"], r["condition"]
    if p not in t.index: continue
    n_kac = t.loc[p, "precursor_KAC" if g == "Precursor" else "normal_KAC"]; n_all = t.loc[p, "precursor_AT2L" if g == "Precursor" else "normal_AT2L"]
    chk.append((n_all - r["n_cells"], n_kac / n_all - r["KAC_frac"]))
chk = np.array(chk); log(f"与第二阶段 KAC 比例表核对：{len(chk)} 个伪批量；细胞数最大差 {np.abs(chk[:,0]).max():.0f}；比例最大差 {np.abs(chk[:,1]).max():.2e}")
usable = t[(t.normal_AT2L > 0) & (t.precursor_AT2L > 0)]
msg = []
msg.append(f"有 Normal 与前体两侧 AT2 谱系细胞的患者 {len(usable)} 名（P7 无前体）")
msg.append(f"KAC：两侧都 >=50 的可配对数 = {int(t.pair_KAC.sum())}；AT2L-other：两侧都 >=50 的可配对数 = {int(t.pair_other.sum())}")
msg.append(f"KAC 细胞数：Normal 侧 中位 {int(usable.normal_KAC.median())}（范围 {int(usable.normal_KAC.min())}–{int(usable.normal_KAC.max())}），前体侧 中位 {int(usable.precursor_KAC.median())}（范围 {int(usable.precursor_KAC.min())}–{int(usable.precursor_KAC.max())}）")
msg.append(f"AT2L-other 细胞数：Normal 侧 中位 {int(usable.normal_other.median())}，前体侧 中位 {int(usable.precursor_other.median())}")
msg.append("KAC 两侧不足 50 的患者：" + ", ".join(f"{p}(Normal {int(r.normal_KAC)}，前体 {int(r.precursor_KAC)})" for p, r in usable[~usable.pair_KAC].iterrows()))
msg.append("AT2L-other 两侧不足 50 的患者：" + (", ".join(f"{p}(Normal {int(r.normal_other)}，前体 {int(r.precursor_other)})" for p, r in usable[~usable.pair_other].iterrows()) or "无"))
for x in msg: log(x)
open(f"{OUT}/state_matched/step0_summary.txt", "w").write("\n".join(msg) + "\n")
pd.set_option("display.width", 250)
log("逐患者：\n" + t[["normal_KAC", "precursor_KAC", "pair_KAC", "normal_other", "precursor_other", "pair_other"]].to_string())
log("===== 第 0 步完成（已停下，等待用户确认）=====")
