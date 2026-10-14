#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_spatial
L=results/mech/rctd; mkdir -p $L
echo "1/3 导出参考……"; python scripts/112a_rctd_prep.py > $L/log_prep.txt 2>&1; tail -2 $L/log_prep.txt
echo "2/3 RCTD 逐切片（56 张，每张几分钟；可中断后重跑续接）……"
caffeinate -i Rscript scripts/112b_rctd.R > $L/log_rctd.txt 2>&1; tail -3 $L/log_rctd.txt
echo "3/3 统计……"; python scripts/113_rctd_stats.py 2>&1 | grep -v Warn | tail -8
