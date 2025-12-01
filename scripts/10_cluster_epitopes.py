import argparse
import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

from Bio.SVDSuperimposer import SVDSuperimposer  # Import for Kabsch algorithm
from Bio import pairwise2  # Import for sequence alignment
import networkx as nx  # Import for graph operations

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# 3-letter to 1-letter amino acid code mapping
AMINO_ACID_3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # Common non-standard residues that might appear, map to 'X'
    "ASX": "X",
    "GLX": "X",
    "UNK": "X",
    "SEC": "U",
    "PYL": "O",
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Cluster epitopes based on structural and sequence similarity."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Path to the input epitope summary CSV file (e.g., data/extract_epitope/epitope_summary.csv).",
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
        help="Directory to save the output clustered CSV file.",
    )
    parser.add_argument(
        "--n_res_threshold",
        type=int,
        default=5,
        help="Maximum allowed difference in residue count for coarse pre-filtering.",
    )
    parser.add_argument(
        "--rmsd_threshold",
        type=float,
        default=2.5,
        help="RMSD threshold (Angstroms) for structural similarity.",
    )
    parser.add_argument(
        "--seq_id_threshold",
        type=float,
        default=40.0,
        help="Sequence Identity threshold (percentage) for sequence similarity.",
    )
    parser.add_argument(
        "--workers", type=int, default=8, help="Number of parallel workers."
    )
    return parser.parse_args()


def load_pocket_data(
    npz_path: Path, input_npz_dir: Path
) -> Optional[Tuple[np.ndarray, str, int]]:
    """
    Loads pocket data from an NPZ file, extracts C-alpha coordinates and sequence.

    Returns:
        Tuple[np.ndarray, str, int]: (C-alpha coordinates, 1-letter amino acid sequence, number of residues).
                                   Returns None if loading fails or no C-alpha atoms found.
    """
    try:
        # If pocket_npz_path is relative, construct absolute path using input_npz_dir
        if not npz_path.is_absolute():
            # Assuming the path in CSV is just the filename relative to input_npz_dir
            full_npz_path = input_npz_dir / npz_path.name
        else:
            full_npz_path = npz_path

        if not full_npz_path.exists():
            logger.warning(f"NPZ file not found: {full_npz_path}")
            return None

        data = np.load(full_npz_path)

        coords = data["coords"]
        atom_names = data["atom_names"]  # ADDED: Load atom names
        res_names = data["res_names"]

        ca_indices = np.where(atom_names == "CA")[0]

        if len(ca_indices) == 0:
            logger.warning(f"No C-alpha atoms found in {full_npz_path}.")
            return None

        residue_indices = data["residue_indices"][ca_indices]

        sort_order = np.argsort(residue_indices)

        ca_coords = coords[ca_indices][sort_order]
        ca_res_names = res_names[ca_indices][sort_order]

        sequence = "".join(
            [AMINO_ACID_3_TO_1.get(res.upper(), "X") for res in ca_res_names]
        )

        num_residues = len(ca_res_names)

        return ca_coords, sequence, num_residues

    except Exception as e:
        logger.error(f"Failed to load or parse NPZ {npz_path}: {e}")
        return None


def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calculates the Root Mean Square Deviation (RMSD) between two sets of coordinates
    using the Kabsch algorithm for optimal superposition.
    Assumes coords1 and coords2 have the same number of points.
    """
    if len(coords1) != len(coords2):
        return float("inf")  # Indicate non-comparable

    sup = SVDSuperimposer()
    sup.set(coords1, coords2)
    sup.run()
    return sup.get_rms()


def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """
    Calculates the sequence identity percentage between two sequences
    using global alignment.
    """
    if not seq1 or not seq2:
        return 0.0

    # Use global alignment with match score 1, mismatch 0, gap open -1, gap extend -0.1
    alignments = pairwise2.align.globalms(seq1, seq2, 1, 0, -1, -0.1)
    if not alignments:
        return 0.0

    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, begin, end = best_alignment

    matches = 0
    # The length of the shorter sequence is typically used as the denominator for identity.
    # Here, we count common non-gapped positions where residues are identical.
    common_length = min(len(seq1), len(seq2))

    if common_length == 0:
        return 0.0

    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i] and aligned_seq1[i] != "-":
            matches += 1

    return (matches / common_length) * 100


def process_pair(
    idx_i: int,
    idx_j: int,
    coords_i: np.ndarray,
    coords_j: np.ndarray,
    seq_i: str,
    seq_j: str,
    rmsd_threshold: float,
    seq_id_threshold: float,
) -> Optional[Tuple[int, int, float, float]]:
    """
    Processes a pair of epitopes, calculates RMSD and Sequence Identity.
    Returns the pair indices, RMSD, and Sequence Identity if thresholds are met.
    """
    rmsd_val = calculate_rmsd(coords_i, coords_j)
    if rmsd_val > rmsd_threshold:
        return None  # RMSD threshold not met

    seq_id_val = calculate_sequence_identity(seq_i, seq_j)
    if seq_id_val < seq_id_threshold:
        return None  # Sequence Identity threshold not met

    return idx_i, idx_j, rmsd_val, seq_id_val


def main():
    args = parse_arguments()

    # Setup directories
    input_csv_path = Path(args.input_csv)
    input_npz_dir = Path(args.input_npz_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting epitope clustering with parameters:")
    logger.info(f"  Input CSV: {input_csv_path}")
    logger.info(f"  Input NPZ Directory: {input_npz_dir}")
    logger.info(f"  Output Directory: {output_dir}")
    logger.info(f"  Residue Count Threshold: {args.n_res_threshold}")
    logger.info(f"  RMSD Threshold: {args.rmsd_threshold}")
    logger.info(f"  Sequence Identity Threshold: {args.seq_id_threshold}")
    logger.info(f"  Workers: {args.workers}")

    # Load Metadata
    if not input_csv_path.exists():
        logger.error(f"Input CSV not found: {input_csv_path}")
        return

    df_epitopes = pd.read_csv(input_csv_path)
    logger.info(f"Loaded {len(df_epitopes)} epitope entries from {input_csv_path}")

    processed_epitopes_data = []

    results_from_loading = Parallel(n_jobs=args.workers)(
        delayed(load_pocket_data)(Path(row["pocket_npz_path"]), input_npz_dir)
        for _, row in tqdm(
            df_epitopes.iterrows(), total=len(df_epitopes), desc="Loading NPZ data"
        )
    )

    for i, result in enumerate(results_from_loading):
        if result:
            ca_coords, sequence, num_residues = result
            original_row = df_epitopes.iloc[i]
            processed_epitopes_data.append(
                {
                    "epitope_id": original_row["epitope_id"],
                    "ca_coords": ca_coords,
                    "sequence": sequence,
                    "num_residues": num_residues,
                    "original_row": original_row,
                }
            )

    if not processed_epitopes_data:
        logger.error("No valid epitope data loaded. Exiting.")
        return

    logger.info(
        f"Successfully loaded data for {len(processed_epitopes_data)} epitopes."
    )

    # 1. Coarse Pre-filtering based on residue count
    valid_pairs = []
    num_epitopes = len(processed_epitopes_data)
    logger.info(
        f"Applying coarse pre-filtering based on residue count (N_res_threshold={args.n_res_threshold})."
    )
    for i in tqdm(range(num_epitopes), desc="Pre-filtering epitope pairs"):
        for j in range(i + 1, num_epitopes):
            res_i = processed_epitopes_data[i]["num_residues"]
            res_j = processed_epitopes_data[j]["num_residues"]

            # If `n_res_threshold` is 0, this acts as a strict length filter for RMSD
            if abs(res_i - res_j) <= args.n_res_threshold:
                valid_pairs.append((i, j))

    logger.info(f"Found {len(valid_pairs)} pairs after coarse pre-filtering.")

    # 2. Structural and Sequence Alignment for valid_pairs
    logger.info("Calculating RMSD and Sequence Identity for valid pairs.")
    similarity_results = Parallel(n_jobs=args.workers)(
        delayed(process_pair)(
            i,
            j,
            processed_epitopes_data[i]["ca_coords"],
            processed_epitopes_data[j]["ca_coords"],
            processed_epitopes_data[i]["sequence"],
            processed_epitopes_data[j]["sequence"],
            args.rmsd_threshold,
            args.seq_id_threshold,
        )
        for i, j in tqdm(valid_pairs, desc="Calculating similarities")
    )

    # Filter out None results (where RMSD or SeqID thresholds are not met)
    filtered_similarities = [res for res in similarity_results if res is not None]
    logger.info(
        f"Found {len(filtered_similarities)} pairs meeting similarity thresholds."
    )

    # 3. Build Similarity Graph
    G = nx.Graph()
    for i in range(num_epitopes):
        G.add_node(i, epitope_id=processed_epitopes_data[i]["epitope_id"])

    for i, j, rmsd_val, seq_id_val in filtered_similarities:
        G.add_edge(i, j, rmsd=rmsd_val, seq_id=seq_id_val)

    logger.info(
        f"Built graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges."
    )

    # 4. Clustering using Connected Components
    clusters = list(nx.connected_components(G))
    logger.info(f"Found {len(clusters)} clusters.")

    # 5. Assign Cluster IDs and Select Representatives
    epitope_cluster_data = []
    cluster_id_counter = 0

    # Initialize all with default cluster_id and is_representative
    for i in range(num_epitopes):
        original_row_dict = processed_epitopes_data[i]["original_row"].to_dict()
        original_row_dict["cluster_id"] = -1  # Not yet clustered or singleton
        original_row_dict["is_representative"] = False
        epitope_cluster_data.append(original_row_dict)

    for cluster in clusters:
        if len(cluster) > 1:  # Only assign cluster_id to non-singleton clusters
            cluster_id_counter += 1
            # Sort cluster members by original_row['pdb'] and then epitope_id for consistent representative selection
            sorted_cluster_members = sorted(
                list(cluster),
                key=lambda idx: (
                    processed_epitopes_data[idx]["original_row"]["pdb_id"],
                    processed_epitopes_data[idx]["epitope_id"],
                ),
            )

            representative_idx = sorted_cluster_members[
                0
            ]  # First member after sorting is the representative

            for member_idx in cluster:
                epitope_cluster_data[member_idx]["cluster_id"] = cluster_id_counter
                if member_idx == representative_idx:
                    epitope_cluster_data[member_idx]["is_representative"] = True
        else:  # Handle singletons. Assign a unique cluster ID too.
            singleton_idx = list(cluster)[0]
            cluster_id_counter += 1
            epitope_cluster_data[singleton_idx]["cluster_id"] = cluster_id_counter
            epitope_cluster_data[singleton_idx]["is_representative"] = (
                True  # Singletons are their own representatives
            )

    # 6. Save Results
    result_df = pd.DataFrame(epitope_cluster_data)
    output_csv_path = output_dir / "epitope_cluster_summary.csv"
    result_df.to_csv(output_csv_path, index=False)

    logger.info(
        f"Completed clustering. Saved {len(result_df)} records to {output_csv_path}"
    )
    logger.info(f"Total clusters (including singletons): {cluster_id_counter}")


if __name__ == "__main__":
    main()
