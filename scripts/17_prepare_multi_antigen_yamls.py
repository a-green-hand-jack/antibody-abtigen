import argparse
import logging
import shutil
import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

import gemmi
import numpy as np
import pandas as pd
import yaml
from Bio import Align

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate multi-antigen YAML files for BoltzGen.")
    parser.add_argument("--input_dir", type=str, default="data/cleaned_split_by_antibody_multi", help="Input directory containing clusters.")
    parser.add_argument("--output_cif_dir", type=str, default="data/antibody_multi", help="Output directory for renamed CIF files.")
    parser.add_argument("--output_yaml_dir", type=str, default="data/yamls/multi_antigen", help="Output directory for generated YAML files.")
    parser.add_argument("--summary_csv", type=str, default="data/raw_data/meta/sabdab_summary_all.tsv", help="Path to SAbDab summary TSV.")
    parser.add_argument("--cdr_csv", type=str, default="data/decuplication/summary-decuplication-distance_threshold_9.csv", help="Path to CDR summary CSV.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility (not actively used in generation but good practice).")
    return parser.parse_args()

def setup_directories(dirs: List[Path]):
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

def load_structure(cif_path: Path) -> Optional[gemmi.Structure]:
    try:
        return gemmi.read_structure(str(cif_path))
    except Exception as e:
        logger.error(f"Failed to read CIF {cif_path}: {e}")
        return None

def get_auth_to_label_mapping(cif_path: Path) -> Dict[str, str]:
    mapping = {}
    try:
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        labels = block.find_loop("_atom_site.label_asym_id")
        auths = block.find_loop("_atom_site.auth_asym_id")
        if labels and auths:
            for l, a in zip(labels, auths):
                if a not in mapping:
                    mapping[a] = l
    except Exception as e:
        logger.error(f"Failed to parse mapping from {cif_path}: {e}")
    return mapping

def extract_ca_coords(structure: gemmi.Structure) -> np.ndarray:
    """
    Extracts CA atom coordinates from all chains in the structure.

    Chains are sorted by length (descending) to ensure consistent ordering
    across different PDB files where chain IDs may differ but chain types
    (Heavy/Light) are the same. Heavy chains are typically longer than Light chains.

    Returns: numpy array of shape (N, 3) where N is the number of CA atoms.
    """
    # First, collect CA coords per chain with chain length
    chain_data = []  # List of (chain_length, ca_coords_list)

    for model in structure:
        for chain in model:
            chain_ca_coords = []
            for res in chain:
                # Only standard amino acids
                try:
                    info = gemmi.find_tabulated_residue(res.name)
                    if not info.is_amino_acid():
                        continue
                except Exception:
                    continue

                # Find CA atom
                ca = res.find_atom("CA", "*")
                if ca:
                    chain_ca_coords.append([ca.pos.x, ca.pos.y, ca.pos.z])

            if chain_ca_coords:
                chain_data.append((len(chain_ca_coords), chain_ca_coords))
        break  # Only first model

    # Sort by chain length (descending) to ensure Heavy chain comes first
    chain_data.sort(key=lambda x: -x[0])

    # Concatenate all CA coords in sorted order
    ca_coords = []
    for _, coords in chain_data:
        ca_coords.extend(coords)

    return np.array(ca_coords) if ca_coords else np.empty((0, 3))


def compute_kabsch_transform(
    source_coords: np.ndarray, target_coords: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the optimal rotation matrix and translation to align source to target.

    Uses the Kabsch algorithm to find R, src_centroid, tgt_centroid such that:
        aligned = (source - src_centroid) @ R + tgt_centroid

    Args:
        source_coords: (N, 3) array of source coordinates
        target_coords: (N, 3) array of target coordinates (must have same N)

    Returns:
        R: (3, 3) rotation matrix
        src_centroid: (3,) source centroid
        tgt_centroid: (3,) target centroid
    """
    assert source_coords.shape == target_coords.shape, "Coordinate arrays must have same shape"
    assert source_coords.shape[0] >= 3, "Need at least 3 points for Kabsch"

    # Compute centroids
    src_centroid = source_coords.mean(axis=0)
    tgt_centroid = target_coords.mean(axis=0)

    # Center the coordinates
    src_centered = source_coords - src_centroid
    tgt_centered = target_coords - tgt_centroid

    # Compute covariance matrix
    H = src_centered.T @ tgt_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation matrix
    R = Vt.T @ U.T

    # Handle reflection case (ensure proper rotation)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    return R, src_centroid, tgt_centroid


def apply_transform_to_structure(
    structure: gemmi.Structure,
    R: np.ndarray,
    src_centroid: np.ndarray,
    tgt_centroid: np.ndarray,
) -> gemmi.Structure:
    """
    Applies Kabsch transform to ALL atoms in the structure.

    Transform: new_pos = (old_pos - src_centroid) @ R + tgt_centroid

    Args:
        structure: gemmi.Structure to transform (modified in place)
        R: (3, 3) rotation matrix
        src_centroid: (3,) source centroid
        tgt_centroid: (3,) target centroid

    Returns:
        The transformed structure (same object, modified in place)
    """
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res:
                    # Get current position
                    pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])

                    # Apply transform
                    new_pos = (pos - src_centroid) @ R + tgt_centroid

                    # Update atom position
                    atom.pos = gemmi.Position(new_pos[0], new_pos[1], new_pos[2])

    return structure


def get_epitope_residues(antibody_path: Path, antigen_path: Path, dist_threshold: float = 4.5) -> Dict[str, List[int]]:
    """
    Identifies antigen residues within distance threshold of antibody.
    Returns: {chain_id (AUTH): [list of 1-based POSITIONAL indices]}
    """
    ab_struc = load_structure(antibody_path)
    ag_struc = load_structure(antigen_path)

    if not ab_struc or not ag_struc:
        return {}

    # Build neighbor search for Antibody
    # +2 buffer for grid efficiency
    ns = gemmi.NeighborSearch(ab_struc[0], ab_struc.cell, dist_threshold + 2).populate()

    epitope_residues = {}

    # Iterate over Antigen atoms
    for model in ag_struc:
        for chain in model:
            chain_id = chain.name # Auth ID
            residues_in_contact = set()

            for i, res in enumerate(chain):
                is_contact = False
                for atom in res:
                    # Check if any antibody atom is close
                    close_atoms = ns.find_neighbors(atom, dist_threshold)
                    if close_atoms:
                        is_contact = True
                        break

                if is_contact:
                    # Use 1-based positional index (i + 1)
                    residues_in_contact.add(i + 1)

            if residues_in_contact:
                epitope_residues[chain_id] = sorted(list(residues_in_contact))

    return epitope_residues
def format_residue_ranges(residues: List[int]) -> str:
    """Converts a list of residue numbers [1, 2, 3, 5] to '1..3,5'."""
    if not residues:
        return ""

    ints = sorted(residues)
    ranges = []
    if not ints:
        return ""

    start = ints[0]
    prev = ints[0]

    for x in ints[1:]:
        if x == prev + 1:
            prev = x
        else:
            if start == prev:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}..{prev}")
            start = x
            prev = x

    if start == prev:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}..{prev}")

    return ",".join(ranges)

def get_masked_ranges(masked_seq: str) -> List[Tuple[int, int]]:
    """
    Finds contiguous ranges of 'X' in the sequence.
    Returns (start_index, end_index) tuples (0-based, inclusive).
    """
    ranges = []
    in_gap = False
    start_gap = -1

    for i, char in enumerate(masked_seq):
        if char == 'X':
            if not in_gap:
                in_gap = True
                start_gap = i
        else:
            if in_gap:
                in_gap = False
                ranges.append((start_gap, i - 1))

    if in_gap:
        ranges.append((start_gap, len(masked_seq) - 1))

    return ranges

def get_structure_sequence(cif_path: Path, chain_id: str) -> Tuple[str, List[int]]:
    """
    Extracts the sequence and corresponding residue IDs from the structure.
    Returns: (sequence_string, list_of_positional_indices_1_based)
    """
    struc = load_structure(cif_path)
    if not struc:
        return "", []

    seq_str = ""
    res_ids = []

    for model in struc:
        for chain in model:
            if chain.name == chain_id:
                for i, res in enumerate(chain):
                    # gemmi.find_tabulated_residue gives standard 1-letter code
                    try:
                        # Only standard residues
                        info = gemmi.find_tabulated_residue(res.name)
                        if info.is_amino_acid(): # or is_nucleic_acid if needed
                            code = info.one_letter_code
                            seq_str += code

                            # Use positional index (1-based)
                            # This ensures alignment maps to the nth residue in the file
                            res_ids.append(i + 1)
                    except:
                        pass
        break

    return seq_str, res_ids

def align_and_map_indices(query_seq: str, target_seq: str, target_indices: List[int], masked_seq: str) -> List[Tuple[int, int]]:
    """
    Aligns query_seq (from CSV) to target_seq (from Structure).
    Maps 'X' regions in masked_seq (same length as query_seq) to target_indices.

    Returns list of (start_residue_id, end_residue_id) inclusive.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local' # Local alignment to handle truncated/extended sequences
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    # Align query (CSV) to target (Structure)
    if not query_seq or not target_seq:
        return []

    alignments = aligner.align(query_seq, target_seq)

    if not alignments:
        return []

    alignment = alignments[0]

    # Iterate through masked ranges in query
    x_ranges = get_masked_ranges(masked_seq)
    mapped_ranges = []

    # We need to map query indices to target INDICES (indices into target_indices list)
    # q_idx -> t_idx
    q_to_t = {}

    # alignment.aligned returns tuple of lists of ranges
    # ([(q_start, q_end), ...], [(t_start, t_end), ...])
    q_ranges, t_ranges = alignment.aligned

    for (q_start, q_end), (t_start, t_end) in zip(q_ranges, t_ranges):
        # Map residues
        # Note: q_end is exclusive
        for i in range(q_end - q_start):
            q_idx = q_start + i
            t_idx = t_start + i
            # target_indices list corresponds to target_seq
            if t_idx < len(target_indices):
                q_to_t[q_idx] = target_indices[t_idx]

    for start, end in x_ranges:

        mapped_in_range = []
        for k in range(start, end + 1):
            if k in q_to_t:
                mapped_in_range.append(q_to_t[k])

        if mapped_in_range:
            p_start = mapped_in_range[0]
            p_end = mapped_in_range[-1]
            mapped_ranges.append((p_start, p_end))

    return mapped_ranges

def main():
    args = parse_arguments()

    input_dir = Path(args.input_dir)
    output_cif_dir = Path(args.output_cif_dir)
    output_yaml_dir = Path(args.output_yaml_dir)

    setup_directories([output_cif_dir, output_yaml_dir])

    logger.info(f"Loading summaries from {args.summary_csv} and {args.cdr_csv}")
    sabdab_df = pd.read_csv(args.summary_csv, sep='\t')
    cdr_df = pd.read_csv(args.cdr_csv)

    clusters = sorted([d for d in input_dir.iterdir() if d.is_dir()])
    logger.info(f"Found {len(clusters)} clusters.")

    for cluster_dir in clusters:
        cluster_id = cluster_dir.name
        samples = [d for d in cluster_dir.iterdir() if d.is_dir()]
        if not samples:
            continue

        pdb_to_sample_dir = {d.name.split('_')[0]: d for d in samples}
        pdbs_in_cluster = list(pdb_to_sample_dir.keys())

        # --- 1. Select Representative Antibody ---
        cluster_meta = sabdab_df[sabdab_df['pdb'].isin(pdbs_in_cluster)].copy()

        if cluster_meta.empty:
            rep_pdb = sorted(pdbs_in_cluster)[0]
            logger.warning(f"Cluster {cluster_id}: No metadata found. Selected {rep_pdb} as fallback.")
        else:
            def get_priority(row):
                organism = str(row['organism']).lower() if pd.notna(row['organism']) else ""
                if 'homo sapiens' in organism: return 0
                if 'mus musculus' in organism: return 1
                return 2

            cluster_meta['priority'] = cluster_meta.apply(get_priority, axis=1)
            cluster_meta = cluster_meta.sort_values(['priority', 'pdb'])
            rep_pdb = cluster_meta.iloc[0]['pdb']
            logger.info(f"Cluster {cluster_id}: Selected {rep_pdb} (Priority {cluster_meta.iloc[0]['priority']})")

        rep_sample_dir = pdb_to_sample_dir.get(rep_pdb)
        if not rep_sample_dir: continue

        # --- 2. Process Representative Antibody (SCAFFOLD) ---
        rep_cdr_rows = cdr_df[cdr_df['pdb'] == rep_pdb]
        if rep_cdr_rows.empty:
            logger.warning(f"Cluster {cluster_id}: No CDR info for {rep_pdb}. Skipping.")
            continue

        cdr_info = rep_cdr_rows.iloc[0]

        # Copy Antibody CIF (Cleaned via Gemmi)
        src_ab_cif = rep_sample_dir / "antibody.cif"
        if not src_ab_cif.exists(): continue

        dst_ab_name = f"{cluster_id}_{rep_pdb}_antibody.cif"
        dst_ab_path = output_cif_dir / dst_ab_name

        # Clean and save
        ab_struc = load_structure(src_ab_cif)
        if ab_struc:
            ab_struc.make_mmcif_document().write_file(str(dst_ab_path))
        else:
            continue

        # Get Auth -> Label mapping
        ab_mapping = get_auth_to_label_mapping(dst_ab_path)

        yaml_entities = []

        chains_to_process = []
        if 'H_chain_id' in cdr_info and pd.notna(cdr_info['H_chain_id']):
            chains_to_process.append(('H', cdr_info['H_chain_id'], cdr_info.get('H_chain_seq'), cdr_info.get('H_chain_masked_seq')))
        if 'L_chain_id' in cdr_info and pd.notna(cdr_info['L_chain_id']):
            chains_to_process.append(('L', cdr_info['L_chain_id'], cdr_info.get('L_chain_seq'), cdr_info.get('L_chain_masked_seq')))

        ab_include = []
        ab_design = []
        ab_structure_groups = []
        ab_exclude = []
        ab_design_insertions = []
        ab_reset = []

        for c_type, c_id_auth, c_seq, c_masked_seq in chains_to_process:
            if not c_id_auth or not c_masked_seq or not isinstance(c_masked_seq, str):
                continue

            # Convert Auth ID to Label ID for YAML
            c_id_label = ab_mapping.get(c_id_auth, c_id_auth)

            ab_include.append({"chain": {"id": c_id_label}})

            # Get Structure Sequence using Auth ID (Gemmi structure uses Auth)
            struc_seq, struc_res_ids = get_structure_sequence(dst_ab_path, c_id_auth)

            # Align and Map CDRs
            query_seq = c_seq if isinstance(c_seq, str) else c_masked_seq.replace('X', '?') # Fallback

            mapped_ranges = align_and_map_indices(query_seq, struc_seq, struc_res_ids, c_masked_seq)

            chain_design_indices = []

            for r_start, r_end in mapped_ranges:
                range_str = f"{r_start}..{r_end}"
                chain_design_indices.append(range_str)
                ab_exclude.append({"chain": {"id": c_id_label, "res_index": range_str}})

                # Insertion length estimation
                orig_len = r_end - r_start + 1
                min_len = max(1, orig_len - 2)
                max_len = orig_len + 5

                ab_design_insertions.append({
                    "insertion": {
                        "id": c_id_label,
                        "res_index": int(r_start),
                        "num_residues": f"{min_len}..{max_len}"
                    }
                })

            if chain_design_indices:
                full_design_str = ",".join(chain_design_indices)
                ab_design.append({"chain": {"id": c_id_label, "res_index": full_design_str}})
                ab_structure_groups.append({"group": {"id": c_id_label, "visibility": 2}})
                ab_structure_groups.append({"group": {"id": c_id_label, "visibility": 0, "res_index": full_design_str}})
                ab_reset.append({"chain": {"id": c_id_label}})
            else:
                ab_structure_groups.append({"group": {"id": c_id_label, "visibility": 2}})

        if ab_include:
            ab_entity = {
                "file": {
                    "path": os.path.relpath(str(dst_ab_path), str(output_yaml_dir)),
                    "include": ab_include,
                    "structure_groups": ab_structure_groups
                }
            }
            if ab_design:
                ab_entity["file"]["design"] = ab_design
            if ab_exclude:
                ab_entity["file"]["exclude"] = ab_exclude
            if ab_design_insertions:
                ab_entity["file"]["design_insertions"] = ab_design_insertions
            if ab_reset:
                ab_entity["file"]["reset_res_index"] = ab_reset

            yaml_entities.append(ab_entity)

        # --- 3. Process Antigens (from ALL samples in cluster) ---
        # First, extract reference antibody CA coordinates for alignment
        ref_ab_struc = load_structure(src_ab_cif)
        if not ref_ab_struc:
            logger.error(f"Cluster {cluster_id}: Failed to load reference antibody for alignment.")
            continue
        ref_ab_coords = extract_ca_coords(ref_ab_struc)

        for pdb, s_dir in pdb_to_sample_dir.items():
            src_ag_cif = s_dir / "antigen.cif"
            if not src_ag_cif.exists():
                continue

            dst_ag_name = f"{cluster_id}_{pdb}_antigen.cif"
            dst_ag_path = output_cif_dir / dst_ag_name

            # Load antigen structure
            ag_struc = load_structure(src_ag_cif)
            if not ag_struc:
                continue

            # For non-representative samples, align antigen to reference coordinate system
            if pdb != rep_pdb:
                # Load this sample's antibody to compute alignment transform
                src_paired_ab = s_dir / "antibody.cif"
                paired_ab_struc = load_structure(src_paired_ab)

                if paired_ab_struc:
                    src_ab_coords = extract_ca_coords(paired_ab_struc)

                    # Check if we have enough matching atoms for alignment
                    if len(src_ab_coords) >= 3 and len(ref_ab_coords) >= 3:
                        # Use minimum of both for alignment
                        min_atoms = min(len(src_ab_coords), len(ref_ab_coords))
                        if len(src_ab_coords) != len(ref_ab_coords):
                            logger.warning(
                                "Cluster %s, PDB %s: Antibody CA count mismatch "
                                "(%d vs %d). Using first %d for alignment.",
                                cluster_id, pdb,
                                len(src_ab_coords), len(ref_ab_coords), min_atoms
                            )

                        # Compute Kabsch transform (source antibody -> reference antibody)
                        R, src_centroid, tgt_centroid = compute_kabsch_transform(
                            src_ab_coords[:min_atoms], ref_ab_coords[:min_atoms]
                        )

                        # Apply transform to the ENTIRE antigen structure
                        # This preserves the antibody-antigen relative position
                        ag_struc = apply_transform_to_structure(
                            ag_struc, R, src_centroid, tgt_centroid
                        )
                        logger.debug(
                            "Cluster %s, PDB %s: Aligned antigen to reference.",
                            cluster_id, pdb
                        )
                    else:
                        logger.warning(
                            "Cluster %s, PDB %s: Not enough CA atoms for alignment "
                            "(src: %d, ref: %d). Antigen will not be aligned.",
                            cluster_id, pdb, len(src_ab_coords), len(ref_ab_coords)
                        )
                else:
                    logger.warning(
                        "Cluster %s, PDB %s: Failed to load paired antibody. "
                        "Antigen will not be aligned.",
                        cluster_id, pdb
                    )

            # Save the (possibly transformed) antigen structure
            ag_struc.make_mmcif_document().write_file(str(dst_ag_path))

            # Get mapping for antigen (after saving, to ensure consistency)
            ag_mapping = get_auth_to_label_mapping(dst_ag_path)

            # Compute epitopes using ORIGINAL (untransformed) structures
            # because epitope detection depends on the original pairing
            src_paired_ab = s_dir / "antibody.cif"
            epitopes = get_epitope_residues(src_paired_ab, src_ag_cif)

            # Reload to get chain info (structure was modified in place)
            ag_struc_for_chains = load_structure(dst_ag_path)
            if not ag_struc_for_chains:
                continue

            # ag_chains are Auth IDs
            ag_chains_auth = [c.name for c in ag_struc_for_chains[0]]
            ag_include = []
            ag_binding = []
            ag_groups = []

            for chain_id_auth in ag_chains_auth:
                chain_id_label = ag_mapping.get(chain_id_auth, chain_id_auth)

                ag_include.append({"chain": {"id": chain_id_label}})
                ag_groups.append({"group": {"id": chain_id_label, "visibility": 2}})

                if chain_id_auth in epitopes:
                    range_str = format_residue_ranges(epitopes[chain_id_auth])
                    if range_str:
                        ag_binding.append({"chain": {"id": chain_id_label, "binding": range_str}})

            ag_entity = {
                "file": {
                    "path": os.path.relpath(str(dst_ag_path), str(output_yaml_dir)),
                    "include": ag_include,
                    "structure_groups": ag_groups,
                }
            }
            if ag_binding:
                ag_entity["file"]["binding_types"] = ag_binding

            yaml_entities.append(ag_entity)

        # --- 4. Write YAML ---
        final_yaml = {"entities": yaml_entities}
        out_yaml_path = output_yaml_dir / f"{cluster_id}.yaml"
        with open(out_yaml_path, 'w') as f:
            yaml.dump(final_yaml, f, sort_keys=False)

        logger.info(f"Generated YAML for cluster {cluster_id} -> {out_yaml_path}")

if __name__ == "__main__":
    main()
