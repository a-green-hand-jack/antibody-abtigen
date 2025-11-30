#!/usr/bin/env python3
"""
Test script for NumpyEpitopeGrouper.

Tests:
1. Load embeddings from HDF5
2. Build index with grouper
3. Compute similarity matrix
4. Find groups with threshold 0.85
5. Save outputs (groups.json, similarity_matrix.h5, grouping_stats.csv)
6. Print report
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import logging
import numpy as np

from antibody_abtigen.epitope_pipeline import (
    HDF5EmbeddingStore,
    NumpyEpitopeGrouper,
    save_groups_json,
    save_grouping_stats_csv,
    generate_grouping_report,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    # Paths
    embeddings_path = project_root / "data/epitope_pipeline/embeddings/test_embeddings.h5"
    output_dir = project_root / "data/epitope_pipeline/grouping"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print("Testing NumpyEpitopeGrouper")
    print(f"{'='*70}\n")

    # Step 1: Load embeddings from HDF5
    print("Step 1: Loading embeddings from HDF5...")
    store = HDF5EmbeddingStore()
    encoder_outputs = store.load_all_encoder_outputs(embeddings_path)

    print(f"  Loaded {len(encoder_outputs)} encoder outputs")
    for output in encoder_outputs[:3]:  # Show first 3
        print(f"    - {output.epitope_id}: {output.total_epitope_residues} epitope residues")
    if len(encoder_outputs) > 3:
        print(f"    ... and {len(encoder_outputs) - 3} more")

    # Step 2: Build index
    print("\nStep 2: Building grouper index...")
    grouper = NumpyEpitopeGrouper(
        similarity_threshold=0.85,
        min_group_size=2,
        exclude_same_pdb=True
    )
    grouper.build_index_from_encoder_outputs(encoder_outputs, use_epitope_embedding=True)
    print(f"  Index built with {len(encoder_outputs)} embeddings")

    # Step 3: Compute similarity matrix
    print("\nStep 3: Computing similarity matrix...")
    similarity_matrix = grouper.compute_similarity_matrix()
    print(f"  Matrix shape: {similarity_matrix.shape}")
    print(f"  Value range: [{similarity_matrix.min():.4f}, {similarity_matrix.max():.4f}]")

    # Show some example similarities
    print(f"\n  Sample pairwise similarities:")
    n = min(5, len(encoder_outputs))
    for i in range(n):
        for j in range(i+1, n):
            sim = similarity_matrix[i, j]
            print(f"    {encoder_outputs[i].epitope_id} vs {encoder_outputs[j].epitope_id}: {sim:.4f}")

    # Step 4: Find groups
    print("\nStep 4: Finding groups with threshold 0.85...")
    grouping_output = grouper.find_groups_detailed(similarity_threshold=0.85)

    print(f"  Total epitopes: {grouping_output.total_epitopes}")
    print(f"  Grouped epitopes: {grouping_output.grouped_epitopes}")
    print(f"  Singletons: {grouping_output.singleton_count}")
    print(f"  Number of groups: {len(grouping_output.groups)}")

    # Step 5: Save outputs
    print("\nStep 5: Saving outputs...")

    # 5a: Save groups JSON
    groups_json_path = output_dir / "groups.json"
    save_groups_json(grouping_output, groups_json_path)
    print(f"  Saved groups.json to {groups_json_path}")

    # 5b: Save similarity matrix (full)
    sim_matrix_path = output_dir / "similarity_matrix.h5"
    grouper.save_similarity_matrix(sim_matrix_path)
    print(f"  Saved similarity_matrix.h5 to {sim_matrix_path}")

    # 5c: Save sparse similarity (only pairs above threshold)
    sparse_sim_path = output_dir / "similarity_sparse.h5"
    grouper.save_sparse_similarity(sparse_sim_path, threshold=0.85)
    print(f"  Saved similarity_sparse.h5 to {sparse_sim_path}")

    # 5d: Save stats CSV
    stats_path = output_dir / "grouping_stats.csv"
    save_grouping_stats_csv(grouping_output, stats_path)
    print(f"  Saved grouping_stats.csv to {stats_path}")

    # Step 6: Print report
    print("\nStep 6: Generating report...")
    report = generate_grouping_report(grouping_output)
    print(report)

    # Save report
    report_path = output_dir / "grouping_report.txt"
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"\nSaved report to {report_path}")

    # Test with lower threshold to see more groups
    print("\n" + "="*70)
    print("Testing with lower threshold (0.70)...")
    print("="*70)

    grouping_output_70 = grouper.find_groups_detailed(similarity_threshold=0.70)
    print(f"  Threshold 0.70: {len(grouping_output_70.groups)} groups, {grouping_output_70.grouped_epitopes} grouped")

    # Save lower threshold results too
    save_groups_json(grouping_output_70, output_dir / "groups_threshold_0.70.json")

    print("\n" + "="*70)
    print("Grouper test completed successfully!")
    print("="*70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
