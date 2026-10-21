"""A6 第 3 步：合并 Python（PHReg，复现原报告）与 R（coxph 诊断、Firth）结果。环境 ggn_target。
用法：python scripts/review_A/A6_assemble.py <A6_prepare 的输出目录>
- A6a/§6.5(a)：每个队列内模型的 n、事件数、参数个数、缺失值处理、收敛、Schoenfeld 全局与逐变量 P；拟合失败列出。
- §6.5(b)：合并分析（纳入队列、各队列估计与 SE、方法、τ²、I²、CI、PI）。
- 核对 Python 与 R 的评分系数一致（同一设计矩阵、同为 Breslow）。"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import check_plan, log, dump
import numpy as np, pandas as pd

sha = check_plan(); log("方案哈希校验通过", sha[:16])
RUN = sys.argv[1].rstrip("/") + "/"
py = pd.read_csv(RUN + "A6_cohort_models_python.tsv", sep="\t")
rc = pd.read_csv(RUN + "R/A6_R_cox.tsv", sep="\t")
z = pd.read_csv(RUN + "R/A6_R_zph_terms.tsv", sep="\t")
pop = pd.read_csv(RUN + "A6_analysis_population.tsv", sep="\t")
M = py.merge(rc, on=["model", "cohort"], how="left", suffixes=("", "_R"))
M["abs_diff_beta_py_vs_R"] = (M["beta"] - M["beta_R"]).abs()
M = M.merge(pop[["cohort", "excluded_missing_or_nonpositive_time", "excluded_no_expression", "age_missing_median_imputed",
                 "sex_unrecognised_coded_as_0"]], on="cohort", how="left")
M.to_csv(RUN + "A6a_cohort_models_full.tsv", sep="\t", index=False)
small = pd.read_csv(RUN + "R/A6_R_small_event.tsv", sep="\t") if os.path.exists(RUN + "R/A6_R_small_event.tsv") else pd.DataFrame()
pooled = pd.read_csv(RUN + "A6_pooled.tsv", sep="\t")
okm = M[M.fit_status == "ok"]
summary = {
    "plan_sha256": sha,
    "n_cohort_level_fits": int(len(M[M.cohort != "ALL"])),
    "n_fit_failed_python": int((M.fit_status.astype(str).str.startswith("failed")).sum()),
    "models_not_reproducible": M.loc[M.cohort == "ALL", ["model", "fit_status"]].astype(str).to_dict("records"),
    "n_R_not_ok": int(((M.cohort != "ALL") & (M.status.astype(str) != "ok")).sum()),
    "R_not_ok_list": M.loc[(M.cohort != "ALL") & (M.status.astype(str) != "ok"), ["model", "cohort", "status", "warnings"]].astype(str).to_dict("records"),
    "max_abs_diff_beta_python_vs_R": float(okm["abs_diff_beta_py_vs_R"].max()) if len(okm) else None,
    "PH_global_p_lt_0.05": M.loc[M.zph_global_p < 0.05, ["model", "cohort", "zph_global_p"]].to_dict("records"),
    "PH_score_term_p_lt_0.05": M.loc[M.zph_score_p < 0.05, ["model", "cohort", "zph_score_p"]].to_dict("records"),
    "small_event_cohorts": small.to_dict("records") if len(small) else [],
    "pooled_not_reproduced": pooled.loc[pooled.reproduced == False, "model"].tolist() if "reproduced" in pooled else None,
    "note": ("PH 诊断仅适用于队列内 Cox；合并估计不做 PH 诊断（方案 6.5）。Schoenfeld 显著只提示该队列该模型可能违背比例风险，"
             "不据此更换模型或删除队列。Firth 与简化模型仅针对事件数 <20 的队列，失败即如实列出。"),
}
dump(summary, RUN + "A6_results.json")
log("A6 汇总完成；PH 全局 P<0.05 的队列内模型数", len(summary["PH_global_p_lt_0.05"]), "；Python 与 R 系数最大差",
    summary["max_abs_diff_beta_python_vs_R"])
