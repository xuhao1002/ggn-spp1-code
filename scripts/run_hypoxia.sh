#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; mkdir -p results/hypoxia
conda activate ggn_target; python scripts/106_hypoxia_H1.py 2>&1 | tee results/hypoxia/log_106.txt | grep -E "^(H1|HIF used|DONE)|Error"
conda activate ggn_spatial; python scripts/107_hypoxia_H2.py 2>&1 | tee results/hypoxia/log_107.txt | grep -E "^(H2|DONE)|Error"
