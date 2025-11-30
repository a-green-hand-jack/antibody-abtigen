#!/usr/bin/env python3
"""
Test script for EpitopePipeline orchestrator.

Tests the full pipeline on a small sample of CIF files.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import logging

from antibody_abtigen.epitope_pipeline import (
    EpitopePipeline,
    create_default_config,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    # Paths
    input_dir = project_root / "data/cleaned_cif"
    output_dir = project_root / "data/epitope_pipeline/orchestrator_test"
    sabdab_path = project_root / "data/meta/sabdab_summary_all.tsv"

    # Use a small sample for testing
    limit = 10

    print(f"\n{'='*70}")
    print("Testing EpitopePipeline Orchestrator")
    print(f"{'='*70}\n")

    print(f"Input: {input_dir}")
    print(f"Output: {output_dir}")
    print(f"Limit: {limit}")
    print()

    # Check input directory
    if not input_dir.exists():
        print(f"Error: Input directory not found: {input_dir}")
        print("Using raw_cif directory instead...")
        input_dir = project_root / "data/raw_cif"
        if not input_dir.exists():
            print(f"Error: raw_cif directory also not found")
            return 1

    # Create config
    config = create_default_config(
        data_dir=output_dir,
        device="cuda",
        similarity_threshold=0.85
    )

    # Initialize pipeline
    pipeline = EpitopePipeline(
        config=config,
        sabdab_summary_path=sabdab_path if sabdab_path.exists() else None,
        verbose=True
    )

    # For testing, we'll skip cleaning if using already cleaned files
    skip_stages = []
    if "cleaned_cif" in str(input_dir):
        skip_stages = ['clean']  # Use pre-cleaned structures

    # Run pipeline
    result = pipeline.run_full(
        input_dir=input_dir,
        output_dir=output_dir,
        limit=limit,
        skip_stages=skip_stages
    )

    # Print result
    print(f"\n{'='*70}")
    print("Pipeline Result")
    print(f"{'='*70}")
    print(f"Success: {result.success}")
    print(f"Stages completed: {result.stages_completed}")
    print(f"Total structures: {result.total_structures}")
    print(f"Cleaned: {result.cleaned_structures}")
    print(f"Embedded: {result.embedded_structures}")
    print(f"Groups: {result.groups_found}")
    print(f"Aligned: {result.groups_aligned}")
    print(f"Time: {result.total_time_seconds:.1f}s")

    if result.errors:
        print(f"\nErrors ({len(result.errors)}):")
        for err in result.errors:
            print(f"  - {err}")

    # Verify output files
    print(f"\nOutput files:")
    expected_files = [
        output_dir / "pipeline_summary.json",
        output_dir / "embeddings" / "embeddings.h5",
        output_dir / "grouping" / "groups.json",
    ]

    for path in expected_files:
        status = "✓" if path.exists() else "✗"
        print(f"  {status} {path.relative_to(output_dir)}")

    # Check aligned groups
    aligned_dir = output_dir / "aligned"
    if aligned_dir.exists():
        groups = list(aligned_dir.glob("group_*"))
        print(f"  ✓ aligned/ ({len(groups)} groups)")

    print(f"\n{'='*70}")
    if result.success:
        print("Orchestrator test PASSED!")
    else:
        print("Orchestrator test FAILED")
    print(f"{'='*70}")

    return 0 if result.success else 1


if __name__ == "__main__":
    sys.exit(main())
