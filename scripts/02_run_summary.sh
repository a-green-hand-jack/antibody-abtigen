#!/bin/bash

# This script runs the summary processing step using the specialized 'anarci_env'.
# It activates the environment and forwards all arguments to the python script.

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
mamba activate /ibex/user/wuj0c/env/mamba/anarci_env
# Run the python script
# Note: We use the python from the activated environment
/ibex/user/wuj0c/env/mamba/anarci_env/bin/python scripts/reference/01_summary.py "$@"
# ./scripts/02_run_summary.sh --summary_dir data/raw_data/meta/sabdab_summary_all.tsv --output_dir data/summary/
