#!/bin/bash
# §12.7 全流程：轮数验证 → 不合并全量 → 合并全量 → 敏感性(10 样本)。中断后重新运行本脚本即可续跑。
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_spatial
export C2L_ACCEL=mps
L=results/30_spatial/logs; mkdir -p $L
nohup caffeinate -i bash -c '
set -e
L=results/30_spatial/logs
C2L_EPOCHS=10000 python scripts/92_c2l_map.py merged 6 GSM9226169_P1_LUAD >> $L/epochcheck_map.txt 2>&1
E=$(python scripts/93_epoch_check.py | tee -a $L/epochcheck.txt | grep USE_EPOCHS | cut -d= -f2)
echo "use epochs $E" >> $L/progress.txt
C2L_EPOCHS=$E python scripts/92_c2l_map.py unmerged 6 all >> $L/unmerged_n6.txt 2>&1;  echo "unmerged done $(date)" >> $L/progress.txt
C2L_EPOCHS=$E python scripts/92_c2l_map.py merged 6 all >> $L/merged_n6.txt 2>&1;      echo "merged done $(date)" >> $L/progress.txt
S="GSM9226172_P3_AIS GSM9226173_P3_LUAD GSM9226178_P5_AIS GSM9226179_P5_LUAD GSM9226180_P6_AAH GSM9226181_P6_LUAD GSM9226184_P8_AIS GSM9226185_P8_LUAD GSM9226189_P10_MIA GSM9226190_P10_LUAD"
for n in 3 10; do C2L_EPOCHS=$E python scripts/92_c2l_map.py merged $n $S >> $L/sens_n$n.txt 2>&1; done; echo "sensitivity done $(date)" >> $L/progress.txt
' > /dev/null 2>&1 &
echo "已在后台启动全流程，进程号 $!"
