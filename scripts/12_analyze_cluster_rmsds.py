import argparse
import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from Bio.SVDSuperimposer import SVDSuperimposer

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze intra-cluster RMSDs.")
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Path to the epitope cluster summary CSV file.",
    )
    parser.add_argument(
        "--input_npz_dir",
        type=str,
        required=True,
        help="Directory containing the epitope pocket NPZ files.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save analysis results (e.g., data/analysis/rmsd).",
    )
    parser.add_argument(
        "--workers", type=int, default=8, help="Number of parallel workers."
    )
    return parser.parse_args()


def load_ca_coords(npz_path: Path) -> Optional[np.ndarray]:
    """
    Loads C-alpha coordinates from an NPZ file.
    """
    try:
        if not npz_path.exists():
            return None
        data = np.load(npz_path)

        coords = data["coords"]
        atom_names = data["atom_names"]
        residue_indices = data["residue_indices"]

        # Filter CA
        ca_mask = atom_names == "CA"
        if not np.any(ca_mask):
            return None

        ca_coords = coords[ca_mask]
        ca_res_indices = residue_indices[ca_mask]

        # Sort by residue index to ensure consistent ordering
        sort_idx = np.argsort(ca_res_indices)
        return ca_coords[sort_idx]

    except Exception as e:
        logger.error(f"Error loading NPZ {npz_path}: {e}")
        return None


def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calculates RMSD between two sets of coordinates using Kabsch algorithm.
    """
    if len(coords1) != len(coords2):
        return float("inf")

    sup = SVDSuperimposer()
    sup.set(coords1, coords2)
    sup.run()
    return sup.get_rms()


def process_cluster_rmsd(
    cluster_id: int,
    cluster_group: pd.DataFrame,
    input_npz_dir: Path,
    output_details_dir: Path,
) -> Optional[Dict]:
    """
    Calculates pairwise RMSDs for a cluster, saves detailed CSV, and returns summary stats.
    """
    members = cluster_group.to_dict("records")
    n_members = len(members)

    if n_members < 2:
        return None

    # 1. Load all coordinates
    # Dictionary mapping epitope_id -> coords
    coords_map = {}
    for member in members:
        epitope_id = member["epitope_id"]
        npz_path = input_npz_dir / Path(member["pocket_npz_path"]).name
        coords = load_ca_coords(npz_path)
        if coords is not None:
            coords_map[epitope_id] = coords
        else:
            logger.warning(
                f"Cluster {cluster_id}: Could not load coords for {epitope_id}"
            )

    valid_epitopes = list(coords_map.keys())
    n_valid = len(valid_epitopes)

    if n_valid < 2:
        return None

    # 2. Pairwise RMSD
    pairwise_results = []
    rmsds = []

    for i in range(n_valid):
        for j in range(i + 1, n_valid):
            epi_i = valid_epitopes[i]
            epi_j = valid_epitopes[j]

            coords_i = coords_map[epi_i]
            coords_j = coords_map[epi_j]

            if len(coords_i) != len(coords_j):
                # Skip incomparable pairs (should represent edge cases if clustering worked well)
                continue

            rmsd_val = calculate_rmsd(coords_i, coords_j)

            pairwise_results.append(
                {"epitope_1": epi_i, "epitope_2": epi_j, "rmsd": rmsd_val}
            )
            rmsds.append(rmsd_val)

    if not pairwise_results:
        return None

    # 3. Save Detailed CSV
    details_df = pd.DataFrame(pairwise_results)
    details_path = output_details_dir / f"cluster_{cluster_id}.csv"
    details_df.to_csv(details_path, index=False)

    # 4. Calculate Stats
    rmsds_arr = np.array(rmsds)
    stats = {
        "cluster_id": cluster_id,
        "num_members": n_members,  # Original members in cluster
        "num_valid_members": n_valid,  # Members with valid coords
        "num_comparisons": len(rmsds),
        "min_rmsd": np.min(rmsds_arr),
        "max_rmsd": np.max(rmsds_arr),
        "mean_rmsd": np.mean(rmsds_arr),
        "median_rmsd": np.median(rmsds_arr),
        "std_rmsd": np.std(rmsds_arr),
    }

    return stats


def main():
    args = parse_arguments()

    input_csv_path = Path(args.input_csv)
    input_npz_dir = Path(args.input_npz_dir)
    output_dir = Path(args.output_dir)
    output_details_dir = output_dir / "details"

    output_dir.mkdir(parents=True, exist_ok=True)
    output_details_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv_path.exists():
        logger.error(f"Input CSV not found: {input_csv_path}")
        return

    df = pd.read_csv(input_csv_path)
    logger.info(f"Loaded {len(df)} records.")

    # Group by cluster and filter singletons
    grouped = df.groupby("cluster_id")
    clusters_to_process = []
    for cluster_id, group in grouped:
        if cluster_id == -1:
            continue
        if len(group) >= 2:
            clusters_to_process.append((cluster_id, group))

    logger.info(f"Found {len(clusters_to_process)} clusters with size >= 2 to analyze.")

    # Parallel Processing
    results = Parallel(n_jobs=args.workers)(
        delayed(process_cluster_rmsd)(
            cluster_id, group, input_npz_dir, output_details_dir
        )
        for cluster_id, group in tqdm(clusters_to_process, desc="Analyzing Clusters")
    )

    # Filter None results
    valid_stats = [r for r in results if r is not None]

    if valid_stats:
        stats_df = pd.DataFrame(valid_stats)
        summary_path = output_dir / "cluster_rmsd_stats.csv"
        stats_df.to_csv(summary_path, index=False)
        logger.info(f"Analysis complete. Summary saved to {summary_path}")
    else:
        logger.warning("No valid cluster stats generated.")


if __name__ == "__main__":
    main()
