"""A2（方案 v1.3 §2-A2）外部单细胞验证的核查、单位表与单位级可视化数据。环境 ggn_target。
A2a 先核查原模型的映射与聚合；A2b 输出单位表与散点数据；核查无误则不新增检验。"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir, dump
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_A2")

# 与 14_microenv_layer1.py 一致的髓系定义
MYE = ["Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil", "Proliferating myeloid"]
CFG = {"GSE164789": dict(unit="patient_histology", tissue_filter="T"),
       "GSE189357": dict(unit="sample", tissue_filter=None)}

report = {"plan_sha256": sha, "cohorts": {}, "discrepancies": []}
for ds, cfg in CFG.items():
    h5 = f"{P}/results/05_validation/{ds}/{ds}_annotated.h5ad"
    a = sc.read_h5ad(h5)                      # 需完整载入以取 SPP1 计数
    log(ds, "载入", a.shape)
    # 与 14_microenv_layer1.py 一致：h5ad 中无 final_label 列，标签由 leiden_1.0 映射注释表得到
    ann = pd.read_csv(f"{P}/results/05_validation/{ds}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
    ob = a.obs.copy()
    lab = ob["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str)
    assert lab.ne("nan").all(), f"{ds}: 存在未注释的 leiden 群"
    # 与 14 一致：一律用原始计数层（.X 为对数归一化值，不能用于 EPCAM<2 规则或伪批量）
    assert "counts" in a.layers, f"{ds}: 缺少 layers['counts']"
    C = a.layers["counts"].tocsr()
    assert "EPCAM" in a.var_names and "SPP1" in a.var_names
    ep = np.asarray(C[:, a.var_names.get_loc("EPCAM")].todense()).ravel()
    keep = lab.isin(MYE).values & (ep < 2)
    if cfg["tissue_filter"]:
        keep &= (ob["tissue"].astype(str) == cfg["tissue_filter"]).values
    log(ds, "髓系且 EPCAM<2", int(keep.sum()), "（tissue 过滤后）")

    ob = ob.loc[keep].copy(); Csub = C[np.where(keep)[0]]
    ob["unit"] = (ob["patient"].astype(str) + "_" + ob["histology"].astype(str)) \
        if cfg["unit"] == "patient_histology" else ob["sample"].astype(str)
    ob["condition"] = np.where(ob["histology"].astype(str) == "IAC", "IAC", "Precursor")
    spp1 = np.asarray(Csub[:, a.var_names.get_loc("SPP1")].todense()).ravel()
    tot = np.asarray(Csub.sum(1)).ravel()                 # 原始计数文库大小
    del a, C

    rows = []
    for u, idx in ob.groupby("unit", observed=True).indices.items():
        n = len(idx)
        rows.append(dict(
            dataset=ds, unit=u,
            patient=str(ob["patient"].iloc[idx[0]]),
            sample=",".join(sorted(ob["sample"].iloc[idx].astype(str).unique())),   # 合并单位列出全部样本
            histology=str(ob["histology"].iloc[idx[0]]),
            condition=str(ob["condition"].iloc[idx[0]]),
            n_cells=int(n),
            included=bool(n >= 50),                      # 原分析的纳入规则
            excluded_reason="" if n >= 50 else "<50 cells",
            # 伪批量尺度（方案 6.1）：先聚合原始计数，再 CP10k，再 log1p
            SPP1_pseudobulk_log1p_cp10k=float(np.log1p(spp1[idx].sum() / max(tot[idx].sum(), 1) * 1e4)),
            SPP1_pos_cell_frac=float((spp1[idx] > 0).mean())))
    U = pd.DataFrame(rows).sort_values(["condition", "unit"])
    U.to_csv(OUT + f"A2_units_{ds}.tsv", sep="\t", index=False)

    inc = U[U.included]
    cnt = inc.condition.value_counts().to_dict()
    report["cohorts"][ds] = {
        "unit_definition": cfg["unit"], "tissue_filter": cfg["tissue_filter"],
        "n_units_total": int(len(U)), "n_units_included": int(len(inc)),
        "n_units_excluded_lt50cells": int((~U.included).sum()),
        "units_by_condition": {k: int(v) for k, v in cnt.items()},
        "n_distinct_patients": int(inc.patient.nunique()),
        "units_per_patient": {int(k): int(v) for k, v in inc.groupby("patient").size().value_counts().items()},
        # 同一患者在两组都有单位：非配对模型把它们当作独立单位，须如实报告（A2a 核查项）
        "patients_with_units_in_both_groups": sorted(
            inc.groupby("patient").condition.nunique().loc[lambda s: s > 1].index.tolist()),
        "design": "pydeseq2 ~ condition（非配对模型，单位视为相互独立）"}
    log(ds, "单位", report["cohorts"][ds]["units_by_condition"],
        "患者", report["cohorts"][ds]["n_distinct_patients"])

# A2a 核查：与稿件已报数字比对（14_console.txt 记录：164789 髓系 15/15；189357 髓系 6/3）
EXPECT = {"GSE164789": {"Precursor": 15, "IAC": 15}, "GSE189357": {"Precursor": 6, "IAC": 3}}
for ds, exp in EXPECT.items():
    got = report["cohorts"][ds]["units_by_condition"]
    if {k: got.get(k, 0) for k in exp} != exp:
        report["discrepancies"].append(
            {"dataset": ds, "expected_units": exp, "observed_units": got,
             "action": "方案 0.3：发现真实差异须查明原因；若原分析有误则重算并更正稿件，原值与更正值同时记录"})
        log("！差异", ds, "期望", exp, "实际", got)
for ds, c in report["cohorts"].items():
    if c["patients_with_units_in_both_groups"]:
        report["discrepancies"].append(
            {"dataset": ds, "type": "same-patient units in both groups",
             "patients": c["patients_with_units_in_both_groups"],
             "note": "稿件称验证队列无同患者配对病灶；非配对 ~ condition 模型将同患者单位视为独立。"
                     "是否构成方案 0.3 所称单位错误、是否重算，由研究者判定后记入偏离记录。"})
        log("！同患者跨组单位", ds, c["patients_with_units_in_both_groups"])
hard = [d for d in report["discrepancies"] if "expected_units" in d]
soft = [d for d in report["discrepancies"] if d.get("type") == "same-patient units in both groups"]
report["conclusion"] = (
    ("单位数与原分析不符：按方案 0.3 查明原因，不得直接沿用原 P 值。" if hard else "单位数与原分析一致。")
    + (" 存在同患者跨组单位（见 discrepancies）：稿件中'验证队列无同患者配对'的表述须更正；"
       "是否另做排除该患者的敏感性分析，由研究者判定并记入偏离记录。" if soft else ""))
dump(report, OUT + "A2_results.json")
log("A2 完成")
