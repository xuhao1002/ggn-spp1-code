#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_target
echo "开始运行 H1c（约 5–10 分钟，期间无输出属正常）……"
python scripts/108_hypoxia_H1c.py > results/hypoxia/log_108.txt 2>&1
echo "运行结束，结果如下："
grep -E "^H1c|DONE|Error|Traceback" results/hypoxia/log_108.txt
tail -3 results/hypoxia/log_108.txt
