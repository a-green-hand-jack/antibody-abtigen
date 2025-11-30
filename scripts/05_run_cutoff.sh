#!/bin/bash

# 这个脚本运行 ./scripts/reference/04_cutoff.py
# 根据时间去除那些在boltz中作为训练数据的蛋白质
# 使用当前工作目录下的 uv 环境

source .venv/bin/activate

uv run python ./scripts/reference/04_cutoff.py --summary_file_path ./data/decuplication/summary-decuplication-distance_threshold_9.csv --cif_dir ./data/raw_data/summary-decuplication-distance_threshold_9 --out_dir ./data/cutoff/ --cutoff_date "2021-09-30"

