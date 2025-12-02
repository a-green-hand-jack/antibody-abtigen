import argparse
import logging
from pathlib import Path
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract clusters with multiple antigens (same antibody)."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        default="data/cluster/antibody_cluster_summary.csv",
        help="Path to the antibody cluster summary CSV file.",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        default="data/cluster/antibody_cluster_multi_antigen.csv",
        help="Path to save the filtered CSV file.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    input_path = Path(args.input_csv)
    output_path = Path(args.output_csv)

    if not input_path.exists():
        logger.error(f"Input CSV not found: {input_path}")
        return

    logger.info(f"Loading {input_path}...")
    df = pd.read_csv(input_path)

    # Calculate cluster sizes
    cluster_counts = df['cluster_id'].value_counts()
    
    # Filter clusters with more than 1 member
    multi_cluster_ids = cluster_counts[cluster_counts > 1].index
    
    logger.info(f"Total clusters: {len(cluster_counts)}")
    logger.info(f"Clusters with > 1 antigen: {len(multi_cluster_ids)}")
    
    df_multi = df[df['cluster_id'].isin(multi_cluster_ids)].copy()
    
    # Sort for better readability: 
    # 1. Cluster ID
    # 2. Representative first (descending order of boolean: True > False)
    # 3. Resolution (best first)
    df_multi.sort_values(
        by=['cluster_id', 'is_representative', 'resolution'], 
        ascending=[True, False, True], 
        inplace=True
    )

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_multi.to_csv(output_path, index=False)
    
    logger.info(f"Saved {len(df_multi)} records to {output_path}")


if __name__ == "__main__":
    main()
