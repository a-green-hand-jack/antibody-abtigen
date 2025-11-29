#!/bin/bash
# =============================================================================
# Quick test script - Run with limited samples before full pipeline
# =============================================================================
#
# Usage:
#   bash scripts/run_test.sh [limit]
#   bash scripts/run_test.sh         # Default: 10 samples
#   bash scripts/run_test.sh 50      # Test with 50 samples
#
# =============================================================================

set -e

LIMIT=${1:-10}
OUTPUT_DIR="./data/multi_antigen/raw"
DATA_DIR="./data/cache"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$DATA_DIR"

echo "=============================================="
echo "Test Run - $LIMIT samples"
echo "=============================================="

# Activate virtual environment
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
else
    echo "ERROR: .venv not found. Run: uv sync"
    exit 1
fi

# Check PyMOL
echo ""
echo "PyMOL status:"
python -c "
from antibody_abtigen import get_pymol_info
info = get_pymol_info()
print(f'  Available: {info[\"available\"]}')
print(f'  Method: {info[\"method\"]}')
if info['python_path']:
    print(f'  Python: {info[\"python_path\"]}')
"

echo ""
echo "Running pipeline with --limit $LIMIT..."
echo ""

antibody-abtigen \
    --output "$OUTPUT_DIR" \
    --data-dir "$DATA_DIR" \
    --limit "$LIMIT"

echo ""
echo "Test complete! Check $OUTPUT_DIR for results."
