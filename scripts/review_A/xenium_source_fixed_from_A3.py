"""重构脚本（2026-10-21）：生成 figures/source_data 所用的 xenium_source_fixed.tsv（表 S3 来源份额、图 4e）。
原表于 2026-10-13 代码勘误（串扰回归量取 min(预期, 上皮观察值)）后交互式生成，脚本未保存。本脚本由 A3 脚本的基础规则输出
（results/review_A3/A3_per_section.tsv，rule=base_ge2、ref_set=all4，其校正已采用 min 规则）推导，并断言与已发布文件一致（容差 1e-6）。"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, check_plan, log, outdir
import pandas as pd, numpy as np
sha = check_plan(); log("方案哈希校验通过", sha[:16])
OUT = outdir("review_xenium_source_fixed")
a = pd.read_csv(f"{P}/results/review_A3/A3_per_section.tsv", sep="\t")
b = a[(a.rule == "base_ge2") & (a.ref_set == "all4")].set_index("section")
fx = pd.read_csv(f"{P}/results/xenium/xenium_per_section.tsv", sep="\t").set_index("section")["fixation"]
t = pd.DataFrame({"group": b.group, "fixation": fx.reindex(b.index), "mac_r": b.mac_share_raw, "epi_r": b.epi_share_raw})
t["oth"] = 1 - t.mac_r - t.epi_r; t["mac_c"] = b.mac_share_corrected; t["epi_c"] = b.epi_share_corrected
t = t.reset_index().rename(columns={"index": "section"})
t.to_csv(OUT + "xenium_source_fixed.tsv", sep="\t", index=False)
ref = pd.read_csv(f"{P}/code_release/source_data/xenium__xenium_source_fixed.tsv", sep="\t").set_index("section")
m = t.set_index("section").loc[ref.index]
assert (m[["group", "fixation"]].astype(str).values == ref[["group", "fixation"]].astype(str).values).all()
d = float(np.abs(m[["mac_r", "epi_r", "oth", "mac_c", "epi_c"]].values - ref[["mac_r", "epi_r", "oth", "mac_c", "epi_c"]].values).max())
assert d < 1e-6, d
log("与已发布文件一致，最大绝对差", d)
