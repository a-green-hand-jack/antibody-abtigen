#!/bin/bash
#SBATCH --job-name=antibody_abtigen
#SBATCH --output=logs/antibody_abtigen_%j.out
#SBATCH --error=logs/antibody_abtigen_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# =============================================================================
# Cross-Species Antibody-Antigen Dataset Builder - IBEX Submission Script
# =============================================================================
#
# Usage:
#   sbatch run_full_pipeline.sh              # Submit to SLURM queue
#   bash run_full_pipeline.sh                # Run interactively (for testing)
#
# Prerequisites:
#   1. Install the package: uv add git+https://github.com/a-green-hand-jack/antibody-abtigen.git
#   2. Create pymol-env: mamba create -n pymol-env -c conda-forge pymol-open-source python=3.11 -y
#   3. Create logs directory: mkdir -p logs
#
# =============================================================================

set -e  # Exit on error

# Configuration
OUTPUT_DIR="./data/multi_antigen/raw"
DATA_DIR="./data/cache"
RESOLUTION_THRESHOLD=2.5
IDENTITY_THRESHOLD=50.0

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$DATA_DIR"
mkdir -p logs

echo "=============================================="
echo "Cross-Species Antibody-Antigen Dataset Builder"
echo "=============================================="
echo ""
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Node: ${SLURM_NODELIST:-$(hostname)}"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  Output directory: $OUTPUT_DIR"
echo "  Data cache: $DATA_DIR"
echo "  Resolution threshold: $RESOLUTION_THRESHOLD Å"
echo "  Identity threshold: $IDENTITY_THRESHOLD%"
echo ""

# Activate virtual environment
if [ -f ".venv/bin/activate" ]; then
    echo "Activating virtual environment..."
    source .venv/bin/activate
else
    echo "ERROR: Virtual environment not found at .venv/bin/activate"
    echo "Please run: uv sync"
    exit 1
fi

# Check if pymol-env exists
echo "Checking PyMOL environment..."
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "WARNING: No conda/mamba found. Will use Biopython for alignment."
    CONDA_CMD=""
fi

if [ -n "$CONDA_CMD" ]; then
    # Check if pymol-env exists
    if $CONDA_CMD env list | grep -q "pymol-env"; then
        echo "Found pymol-env conda environment"
        USE_PYMOL=""
    else
        echo "WARNING: pymol-env not found. Creating it now..."
        $CONDA_CMD create -n pymol-env -c conda-forge pymol-open-source python=3.11 -y
        echo "pymol-env created successfully"
        USE_PYMOL=""
    fi
else
    USE_PYMOL="--no-pymol"
fi

# Verify PyMOL detection
echo ""
echo "Verifying PyMOL setup..."
python -c "
from antibody_abtigen import get_pymol_info, setup_pymol
info = get_pymol_info()
print(f'PyMOL status: {info}')
if not info['available']:
    success, msg = setup_pymol(auto_create=False)
    print(f'Setup result: {msg}')
"

# Run the pipeline
echo ""
echo "=============================================="
echo "Starting pipeline..."
echo "=============================================="
echo ""

antibody-abtigen \
    --output "$OUTPUT_DIR" \
    --data-dir "$DATA_DIR" \
    --resolution "$RESOLUTION_THRESHOLD" \
    --identity "$IDENTITY_THRESHOLD" \
    $USE_PYMOL

# Print summary
echo ""
echo "=============================================="
echo "Pipeline Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  Summary CSV: $OUTPUT_DIR/dataset_summary.csv"
echo "  Processing log: $OUTPUT_DIR/processing_log.json"
echo ""

# Count results
if [ -f "$OUTPUT_DIR/dataset_summary.csv" ]; then
    TOTAL=$(wc -l < "$OUTPUT_DIR/dataset_summary.csv")
    TOTAL=$((TOTAL - 1))  # Subtract header
    SUCCESS=$(grep -c ",success," "$OUTPUT_DIR/dataset_summary.csv" || echo "0")
    FAILED=$(grep -c ",failed," "$OUTPUT_DIR/dataset_summary.csv" || echo "0")

    echo "Results:"
    echo "  Total processed: $TOTAL"
    echo "  Successful: $SUCCESS"
    echo "  Failed: $FAILED"
    echo ""
    echo "Data point folders: $(ls -d "$OUTPUT_DIR"/DP_* 2>/dev/null | wc -l)"
fi
