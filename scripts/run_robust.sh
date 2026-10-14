#!/bin/bash
cd ~/全部数据库/GGN_target
eval "$(conda shell.bash hook)"; conda activate ggn_target
mkdir -p results/sens_protein results/xenium
python scripts/103_xenium_spillover.py 2>&1 | tee results/xenium/log_103.txt | tail -32
python scripts/102_prognosis_robustness.py 2>&1 | tee results/sens_protein/log_102.txt | grep -E "^(main|R[123]|DONE)|Error"
