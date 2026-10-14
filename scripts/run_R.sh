#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_target
python -c "import decoupler" 2>/dev/null || pip install -q decoupler omnipath
echo "开始运行 TF 活性分析（约 10–20 分钟，期间无输出属正常）……"
python scripts/110_tf_activity.py > results/mech/log_110.txt 2>&1
echo "运行结束，结果如下："
grep -vE "^\s*$|Warning|warn" results/mech/log_110.txt | tail -25
