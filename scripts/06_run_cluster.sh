#!/bin/bash

# This script runs the clustering step using the specialized 'mmseqs' environment.
# It activates the environment and runs the python script with specified arguments.

# Source conda setup (adjust path if conda is installed elsewhere)
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "Error: Could not find conda profile script."
    exit 1
fi

eval "$(mamba shell hook --shell bash)"

# Activate the environment
mamba activate /ibex/user/wuj0c/envs/mmseqs

# Run the python script
# We use the python from the activated environment
# We explicitly set the mmseqs executable to 'mmseqs' assuming it's in the path of the activated environment
/ibex/user/wuj0c/envs/mmseqs/bin/python scripts/reference/05_cluster.py \
    --processed_csv_path data/decuplication/summary-decuplication-distance_threshold_9.csv \
    --out_dir data/cluster/ \
    --mmseqs_path mmseqs \
    "$@"
