import argparse
import logging
import shutil
from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Copy cleaned data for strictly multi-antigen clusters to a separate directory."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        default="data/cluster/antibody_cluster_multi_antigen.csv",
        help="Path to the filtered multi-antigen cluster CSV (used for candidate cluster IDs).",
    )
    parser.add_argument(
        "--source_dir",
        type=str,
        default="data/cleaned_split_by_antibody",
        help="Source directory containing all cleaned data.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/cleaned_split_by_antibody_multi",
        help="Target directory to copy the multi-antigen data to.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    input_csv = Path(args.input_csv)
    source_base = Path(args.source_dir)
    output_base = Path(args.output_dir)

    if not input_csv.exists():
        logger.error(f"Input CSV not found: {input_csv}")
        return
    if not source_base.exists():
        logger.error(f"Source directory not found: {source_base}")
        return

    # 1. Clean output directory first to avoid stale data
    if output_base.exists():
        logger.info(f"Cleaning existing output directory: {output_base}")
        shutil.rmtree(output_base)
    output_base.mkdir(parents=True, exist_ok=True)

    # 2. Get Candidate Cluster IDs
    logger.info(f"Loading candidates from {input_csv}...")
    df = pd.read_csv(input_csv)
    candidate_cluster_ids = df['cluster_id'].unique()
    
    logger.info(f"Checking {len(candidate_cluster_ids)} candidate clusters...")

    copied_clusters = 0
    total_epitopes_copied = 0

    # 3. Filter and Copy
    for cluster_id in tqdm(candidate_cluster_ids, desc="Filtering and Copying"):
        cid_str = str(cluster_id)
        src_cluster_dir = source_base / cid_str
        
        if not src_cluster_dir.exists():
            continue
            
        # Count valid subdirectories (epitopes) in the source cluster folder
        # Filter out files like .DS_Store if any, though iterdir usually returns mix
        epitope_dirs = [p for p in src_cluster_dir.iterdir() if p.is_dir()]
        
        if len(epitope_dirs) > 1:
            # STRICT CHECK: Only copy if actual cleaned count > 1
            dst_cluster_dir = output_base / cid_str
            
            # Copy the whole cluster directory
            # shutil.copytree creates the directory
            try:
                shutil.copytree(src_cluster_dir, dst_cluster_dir)
                copied_clusters += 1
                total_epitopes_copied += len(epitope_dirs)
            except Exception as e:
                logger.error(f"Failed to copy cluster {cid_str}: {e}")

    logger.info("-" * 30)
    logger.info(f"Extraction complete.")
    logger.info(f"Clusters processed: {len(candidate_cluster_ids)}")
    logger.info(f"Multi-antigen clusters found (post-cleaning): {copied_clusters}")
    logger.info(f"Total epitopes copied: {total_epitopes_copied}")
    logger.info(f"Output directory: {output_base}")


if __name__ == "__main__":
    main()
