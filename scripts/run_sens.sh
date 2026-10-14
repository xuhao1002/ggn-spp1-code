#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_target
python scripts/98_niche_sensitivity.py 2>&1 | tee results/sens_protein/log_98.txt | grep -v -E '^\s*$' | tail -15
python scripts/99_il1b_power.py 2>&1 | tee results/sens_protein/log_99.txt | tail -14
