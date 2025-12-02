import argparse
import logging
from pathlib import Path
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Cluster epitopes by identical antibody sequences (VH + VL)."
    )
    parser.add_argument(
        "--input_epitope_csv",
        type=str,
        required=True,
        help="Path to the epitope summary CSV file (e.g., data/extract_epitope/epitope_summary.csv).",
    )
    parser.add_argument(
        "--input_summary_csv",
        type=str,
        required=True,
        help="Path to the SAbDab summary cache CSV file containing sequences (e.g., data/summary/summary.csv).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the cluster summary CSV file (e.g., data/cluster).",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    epitope_csv_path = Path(args.input_epitope_csv)
    summary_csv_path = Path(args.input_summary_csv)
    output_dir = Path(args.output_dir)

    if not epitope_csv_path.exists():
        logger.error(f"Epitope CSV not found: {epitope_csv_path}")
        return
    if not summary_csv_path.exists():
        logger.error(f"Summary CSV not found: {summary_csv_path}")
        return

    # Load Data
    logger.info("Loading datasets...")
    df_epitope = pd.read_csv(epitope_csv_path)
    df_summary = pd.read_csv(summary_csv_path)

    logger.info(f"Loaded {len(df_epitope)} epitope records.")
    logger.info(f"Loaded {len(df_summary)} summary records.")

    # Prepare Summary for Merge
    # We need to map pdb, H, L to sequences.
    # summary.csv columns: pdb, H_chain_id, L_chain_id, H_chain_seq, L_chain_seq, resolution
    
    # Select relevant columns to avoid conflicts
    df_seqs = df_summary[[
        'pdb', 'H_chain_id', 'L_chain_id', 
        'H_chain_seq', 'L_chain_seq', 'resolution'
    ]].copy()

    # Handle duplicates in summary if any (same pdb/H/L should be unique usually, but just in case)
    df_seqs.drop_duplicates(subset=['pdb', 'H_chain_id', 'L_chain_id'], inplace=True)

    # Merge
    # epitope: pdb_id, H_chain, L_chain
    # summary: pdb, H_chain_id, L_chain_id
    logger.info("Merging epitope data with antibody sequences...")
    
    merged_df = pd.merge(
        df_epitope,
        df_seqs,
        left_on=['pdb_id', 'H_chain', 'L_chain'],
        right_on=['pdb', 'H_chain_id', 'L_chain_id'],
        how='inner'
    )

    if len(merged_df) < len(df_epitope):
        logger.warning(
            f"Dropped {len(df_epitope) - len(merged_df)} records due to missing sequences in summary.csv."
        )

    # Cluster by Sequences
    logger.info("Clustering by identical VH and VL sequences...")
    
    # Create a unique signature for grouping
    # Handle NaNs in sequences if any (shouldn't be after inner join but good practice)
    merged_df['H_chain_seq'] = merged_df['H_chain_seq'].fillna('')
    merged_df['L_chain_seq'] = merged_df['L_chain_seq'].fillna('')
    
    # Group
    # usage of 'ngroup' provides a dense cluster ID
    merged_df['cluster_id'] = merged_df.groupby(['H_chain_seq', 'L_chain_seq']).ngroup()

    num_clusters = merged_df['cluster_id'].nunique()
    logger.info(f"Identified {num_clusters} unique antibody clusters.")

    # Determine Representative
    # Strategy: Sort by Resolution (asc), then by PDB ID (to be deterministic)
    # Lower resolution value is better (higher quality).
    
    merged_df.sort_values(by=['cluster_id', 'resolution', 'pdb_id'], ascending=[True, True, True], inplace=True)
    
    # Mark representative
    merged_df['is_representative'] = False
    # The first entry in each cluster (after sorting) is the representative
    merged_df.loc[merged_df.groupby('cluster_id').head(1).index, 'is_representative'] = True

    # Prepare Output
    # We need to ensure columns required by downstream scripts are present.
    # scripts/11_align_structures.py needs: cluster_id, epitope_id, pdb_path, pocket_npz_path, is_representative
    # scripts/13_clean_and_split.py needs: cluster_id, epitope_id, H_chain, L_chain, antigen_chains
    
    output_cols = [
        'cluster_id', 'epitope_id', 'pdb_id', 
        'H_chain', 'L_chain', 'antigen_chains',
        'pdb_path', 'pocket_npz_path', 'is_representative',
        'resolution', 'H_chain_seq', 'L_chain_seq' # Extra info
    ]
    
    final_df = merged_df[output_cols]

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "antibody_cluster_summary.csv"
    final_df.to_csv(output_path, index=False)
    
    logger.info(f"Saved cluster summary to {output_path}")
    
    # Print stats
    cluster_sizes = final_df['cluster_id'].value_counts()
    logger.info(f"Cluster size stats: Mean={cluster_sizes.mean():.2f}, Max={cluster_sizes.max()}, Min={cluster_sizes.min()}")


if __name__ == "__main__":
    main()
