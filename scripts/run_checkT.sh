#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_target
python scripts/104_checkT_xenium.py 2>&1 | tee results/xenium/log_104.txt | tail -60
python scripts/105_checkT_hypoxia.py 2>&1 | tee results/sens_protein/log_105.txt | tail -60
