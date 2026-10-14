#!/bin/bash
# 反卷积：可区分性模拟 + 单样本计时（GPU）。用法：bash scripts/run_deconv_bench.sh
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_spatial
export C2L_ACCEL=mps
nohup caffeinate -i bash -c '
python scripts/91_identifiability_sim.py > results/30_spatial/log_91_sim.txt 2>&1 &&
python scripts/92_c2l_map.py merged 6 GSM9226169_P1_LUAD > results/30_spatial/log_92_benchmark.txt 2>&1
' > /dev/null 2>&1 &
echo "已在后台启动，进程号 $!"
