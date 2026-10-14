#!/bin/bash
# 补充数据集。中断或有未完成文件时重新运行即可续传。
# 用法：bash ~/全部数据库/GGN_target/scripts/download_geo_extra.sh
BASE="$HOME/全部数据库/GGN_target/raw"
source "$(dirname "$0")/geo_lib.sh"
get_series GSE164789 scRNA_GSE164789_precursor      # 前体病变单细胞 + scTCR，78 个样本
get_series GSE319666 bulk_GSE319666_HTAN_PML        # HTAN 前体病变 bulk RNA-seq，177 个样本
get_series GSE300827 ST_GSE300827_AIS_MIA           # AIS/MIA 空间转录组，20 个样本，约 678 MB
get_series GSE260770 bulk_GSE260770_sintilimab_GGO  # 信迪利单抗治疗高危磨玻璃结节试验，50 个样本
get_series GSE288479 multi_GSE288479_partsolid      # 部分实性结节磨玻璃区与实性区，38 个样本
get_series GSE203360 scRNA_GSE203360_mGGO           # 混合磨玻璃结节两区配对单细胞，6 个样本
get_series GSE282617 bulk_GSE282617_AIS_MIA_IAC     # AIS→MIA→IAC 进展队列，70 个样本
get_series GSE169033 nanostring_GSE169033           # 免疫基因 panel，85 个样本
get_series GSE52248  bulk_GSE52248_AIS_FFPE         # 正常/AIS/浸润癌 FFPE，18 个样本
echo "== 补充数据处理完毕"
