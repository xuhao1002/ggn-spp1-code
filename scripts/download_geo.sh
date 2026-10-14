#!/bin/bash
# Peng 2026 主数据：单核转录组 + 空间转录组。中断或有未完成文件时重新运行即可续传。
# 用法：bash ~/全部数据库/GGN_target/scripts/download_geo.sh
BASE="$HOME/全部数据库/GGN_target/raw"
source "$(dirname "$0")/geo_lib.sh"
get_series GSE308103 snRNA_GSE308103   # 单核转录组，约 1.56 GB
get_series GSE307534 ST_GSE307534      # 空间转录组，约 10.1 GB
