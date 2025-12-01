import argparse
import logging
from pathlib import Path
from typing import Optional

import gemmi
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
    parser = argparse.ArgumentParser(
        description="Align structures within clusters based on epitope superposition."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Path to the epitope cluster summary CSV file (e.g., data/cluster/epitope_cluster_summary.csv).",
    )
    parser.add_argument(
        "--input_npz_dir",
        type=str,
        required=True,
        help="Directory containing the epitope pocket NPZ files (e.g., data/extract_epitope/pockets).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Base directory to save aligned CIF files (e.g., data/aligned).",
    )
    parser.add_argument(
        "--min_size",
        type=int,
        default=2,
        help="Minimum cluster size to process. Clusters smaller than this will be skipped (or just copied if 1).",
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


def get_transform_from_coords(
    fixed_coords: np.ndarray, mobile_coords: np.ndarray
) -> gemmi.Transform:
    """
    Calculates the transformation (rotation & translation) to align mobile_coords to fixed_coords.
    Returns a gemmi.Transform object.
    """
    sup = SVDSuperimposer()
    sup.set(fixed_coords, mobile_coords)
    sup.run()

    rot, tran = sup.get_rotran()

    # Convert to gemmi.Transform
    # Gemmi Transform matrix is a list of lists [[r11, r12, r13], ...]
    # Gemmi Transform vector is a gemmi.Vec3

    transform = gemmi.Transform()
    transform.mat.fromlist(rot.tolist())
    transform.vec.fromlist(tran.tolist())

    return transform


def process_cluster(
    cluster_id: int,
    cluster_group: pd.DataFrame,
    input_npz_dir: Path,
    output_base_dir: Path,
):
    """
    Aligns all members of a cluster to the representative member.
    """
    # Create cluster directory
    cluster_dir = output_base_dir / str(cluster_id)
    cluster_dir.mkdir(parents=True, exist_ok=True)

    # Identify Representative
    rep_rows = cluster_group[cluster_group["is_representative"]]
    if rep_rows.empty:
        # Fallback: use the first one if no representative marked (shouldn't happen with correct upstream)
        rep_row = cluster_group.iloc[0]
    else:
        rep_row = rep_rows.iloc[0]

    # Load Representative Data
    rep_npz_path = input_npz_dir / Path(rep_row["pocket_npz_path"]).name
    rep_ca_coords = load_ca_coords(rep_npz_path)

    if rep_ca_coords is None:
        logger.warning(
            f"Cluster {cluster_id}: Failed to load representative coordinates. Skipping cluster."
        )
        return

    # Read Representative CIF (full structure)
    try:
        rep_cif_path = Path(rep_row["pdb_path"])
        rep_doc = gemmi.cif.read_file(str(rep_cif_path))
        # Save Representative "As Is"
        out_path = cluster_dir / f"{rep_row['epitope_id']}.cif"
        rep_doc.write_file(str(out_path))
    except Exception as e:
        logger.warning(
            f"Cluster {cluster_id}: Failed to read/write representative CIF {rep_row['pdb_path']}: {e}"
        )
        return

    # Process all members
    for _, row in cluster_group.iterrows():
        # Skip if it is the representative (already handled)
        if row["epitope_id"] == rep_row["epitope_id"]:
            continue

        member_npz_path = input_npz_dir / Path(row["pocket_npz_path"]).name
        member_ca_coords = load_ca_coords(member_npz_path)

        if member_ca_coords is None:
            logger.warning(
                f"Cluster {cluster_id}: Failed to load coords for {row['epitope_id']}. Skipping."
            )
            continue

        # Check if atoms match
        if len(member_ca_coords) != len(rep_ca_coords):
            logger.warning(
                f"Cluster {cluster_id}: Atom count mismatch ({row['epitope_id']}: {len(member_ca_coords)} vs Rep: {len(rep_ca_coords)}). Skipping."
            )
            continue

        # Calculate Transform
        try:
            transform = get_transform_from_coords(rep_ca_coords, member_ca_coords)
        except Exception as e:
            logger.error(
                f"Cluster {cluster_id}: alignment failed for {row['epitope_id']}: {e}"
            )
            continue

        # Apply Transform to Full CIF
        try:
            cif_path = Path(row['pdb_path'])
            # Use read_structure to get a high-level Structure object which supports transform on Models
            st = gemmi.read_structure(str(cif_path))
            
            for model in st:
                model.transform_pos_and_adp(transform)
                
            out_path = cluster_dir / f"{row['epitope_id']}.cif"
            st.make_mmcif_document().write_file(str(out_path))
            
        except Exception as e:
            logger.warning(f"Cluster {cluster_id}: Failed to transform/save {row['epitope_id']}: {e}")
            continue


def main():
    args = parse_arguments()

    input_csv_path = Path(args.input_csv)
    input_npz_dir = Path(args.input_npz_dir)
    output_dir = Path(args.output_dir)

    if not input_csv_path.exists():
        logger.error(f"Input CSV not found: {input_csv_path}")
        return

    df = pd.read_csv(input_csv_path)
    logger.info(f"Loaded {len(df)} records.")

    # Group by cluster
    grouped = df.groupby("cluster_id")

    # Filter clusters by size
    clusters_to_process = []
    for cluster_id, group in grouped:
        if cluster_id == -1:  # Skip unclustered if any
            continue
        if len(group) >= args.min_size:
            clusters_to_process.append((cluster_id, group))

    logger.info(
        f"Found {len(clusters_to_process)} clusters with size >= {args.min_size} to align."
    )

    # Parallel Processing
    Parallel(n_jobs=args.workers)(
        delayed(process_cluster)(cluster_id, group, input_npz_dir, output_dir)
        for cluster_id, group in tqdm(clusters_to_process, desc="Aligning Clusters")
    )

    logger.info("Alignment complete.")


if __name__ == "__main__":
    main()
