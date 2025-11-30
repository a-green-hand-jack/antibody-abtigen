#!/bin/bash

# This script runs the split processing step using the project's uv environment.
# It runs the python script with specified arguments.

# Use uv to run the script, which automatically handles the environment and dependencies
uv run scripts/reference/06_split.py \
    --clu_result_path data/cluster/decuplicated_threshold_9/cluster_results_seq_identity_threshold_0.5/heavy_chain_cdr3/cluster_heavy_chain_cdr3_cluster.tsv \
    --pdb_before_cutoff_in_sabdab_path data/cutoff/before_cutoff_in_sabdab.json \
    --precessed_data_path data/decuplication/summary-decuplication-distance_threshold_9.csv \
    --output_dir ./data/split/ \
    "$@"
