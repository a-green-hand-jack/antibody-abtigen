import argparse
import json
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import gemmi
import pandas as pd
import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate BoltzGen-like YAML files for test entries from raw PDB files.")
    parser.add_argument("--test_entries", type=str, default="data/split/test_entry.json", help="Path to test entry JSON file.")
    parser.add_argument("--pdb_dir", type=str, default="data/raw_data/chothia", help="Directory containing raw PDB files.")
    parser.add_argument("--cdr_csv", type=str, default="data/decuplication/summary-decuplication-distance_threshold_9.csv", help="Path to CDR summary CSV.")
    parser.add_argument("--output_yaml_dir", type=str, default="data/yamls/boltzgen_like/test_entries", help="Output directory for generated YAML files.")
    return parser.parse_args()


def setup_directories(dirs: List[Path]):
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def load_test_entries(test_entries_path: Path) -> List[str]:
    """Load test entry identifiers from JSON file."""
    try:
        with open(test_entries_path, 'r') as f:
            entries = json.load(f)
        return entries
    except Exception as e:
        logger.error(f"Failed to load test entries from {test_entries_path}: {e}")
        return []


def load_structure(pdb_path: Path) -> Optional[gemmi.Structure]:
    try:
        return gemmi.read_structure(str(pdb_path))
    except Exception as e:
        logger.error(f"Failed to read PDB {pdb_path}: {e}")
        return None


def get_auth_to_label_mapping(pdb_path: Path) -> Dict[str, str]:
    """Maps Auth chain IDs to Label chain IDs from PDB file."""
    mapping = {}
    try:
        doc = gemmi.cif.read_file(str(pdb_path))
        block = doc.sole_block()
        labels = block.find_loop("_atom_site.label_asym_id")
        auths = block.find_loop("_atom_site.auth_asym_id")
        if labels and auths:
            for l, a in zip(labels, auths):
                if a not in mapping:
                    mapping[a] = l
    except Exception as e:
        # PDB files may not have the CIF structure, use simple mapping
        pass
    return mapping


def get_epitope_residues(antibody_path: Path, antigen_path: Path, ab_chain: str, ag_chain: str, dist_threshold: float = 4.5) -> List[int]:
    """
    Identifies antigen residues within distance threshold of antibody.
    Returns: [list of 1-based POSITIONAL indices]
    """
    ab_struc = load_structure(antibody_path)
    ag_struc = load_structure(antigen_path)

    if not ab_struc or not ag_struc:
        return []

    # Build neighbor search for Antibody
    ns = gemmi.NeighborSearch(ab_struc[0], ab_struc.cell, dist_threshold + 2).populate()

    residues_in_contact = set()

    # Iterate over Antigen atoms
    for model in ag_struc:
        for chain in model:
            if chain.name != ag_chain:
                continue

            for i, res in enumerate(chain):
                is_contact = False
                for atom in res:
                    close_atoms = ns.find_neighbors(atom, dist_threshold)
                    if close_atoms:
                        is_contact = True
                        break

                if is_contact:
                    residues_in_contact.add(i + 1)

    return sorted(list(residues_in_contact))


def format_residue_ranges(residues: List[int]) -> str:
    """Converts a list of residue numbers [1, 2, 3, 5] to '1..3,5'."""
    if not residues:
        return ""

    ints = sorted(residues)
    ranges = []

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


def main():
    args = parse_arguments()

    pdb_dir = Path(args.pdb_dir)
    output_yaml_dir = Path(args.output_yaml_dir)

    setup_directories([output_yaml_dir])

    # Load test entries
    logger.info(f"Loading test entries from {args.test_entries}")
    test_entries = load_test_entries(Path(args.test_entries))
    logger.info(f"Found {len(test_entries)} test entries")

    # Load CDR information
    logger.info(f"Loading CDR information from {args.cdr_csv}")
    cdr_df = pd.read_csv(args.cdr_csv)

    processed_count = 0
    skipped_count = 0

    for entry in test_entries:
        # Parse entry format: e.g., "9ixv_H_L_A" or "9j3j_B__A"
        # Format: PDB_code_H_chain_L_chain_antigen_chains
        parts = entry.split('_')
        pdb_code = parts[0]
        h_chain = parts[1] if len(parts) > 1 and parts[1] else None
        l_chain = parts[2] if len(parts) > 2 and parts[2] else None
        antigen_chains = parts[3] if len(parts) > 3 else None

        # Load PDB file
        pdb_path = pdb_dir / f"{pdb_code}.pdb"
        if not pdb_path.exists():
            logger.warning(f"Skipping {entry}: PDB file not found at {pdb_path}")
            skipped_count += 1
            continue

        # Load structure
        structure = load_structure(pdb_path)
        if not structure:
            logger.warning(f"Skipping {entry}: Failed to load PDB structure")
            skipped_count += 1
            continue

        # Get mapping
        mapping = get_auth_to_label_mapping(pdb_path)

        # Try to get CDR info from CSV
        cdr_rows = cdr_df[cdr_df['pdb'] == pdb_code]
        if cdr_rows.empty:
            logger.warning(f"Skipping {entry}: No CDR information found for {pdb_code}")
            skipped_count += 1
            continue

        cdr_info = cdr_rows.iloc[0]

        yaml_entities = []

        # --- Process Antibody ---
        ab_include = []
        ab_design = []
        ab_structure_groups = []
        ab_exclude = []
        ab_design_insertions = []
        ab_reset = []

        chains_to_process = []

        # Use chains from test entry if specified
        if h_chain and h_chain != '':
            chains_to_process.append(('H', h_chain, cdr_info.get('H_chain_seq'), cdr_info.get('H_chain_masked_seq')))
        if l_chain and l_chain != '':
            chains_to_process.append(('L', l_chain, cdr_info.get('L_chain_seq'), cdr_info.get('L_chain_masked_seq')))

        # Fallback to CSV if not specified
        if not chains_to_process:
            if pd.notna(cdr_info.get('H_chain_id')):
                chains_to_process.append(('H', cdr_info['H_chain_id'], cdr_info.get('H_chain_seq'), cdr_info.get('H_chain_masked_seq')))
            if pd.notna(cdr_info.get('L_chain_id')):
                chains_to_process.append(('L', cdr_info['L_chain_id'], cdr_info.get('L_chain_seq'), cdr_info.get('L_chain_masked_seq')))

        for c_type, c_id_auth, c_seq, c_masked_seq in chains_to_process:
            if not c_id_auth or not c_masked_seq or not isinstance(c_masked_seq, str):
                continue

            # Convert Auth ID to Label ID
            c_id_label = mapping.get(c_id_auth, c_id_auth)
            ab_include.append({"chain": {"id": c_id_label}})

            # Parse masked sequence to find 'X' regions
            x_ranges = []
            in_gap = False
            start_gap = -1

            for i, char in enumerate(c_masked_seq):
                if char == 'X':
                    if not in_gap:
                        in_gap = True
                        start_gap = i
                else:
                    if in_gap:
                        in_gap = False
                        x_ranges.append((start_gap, i - 1))

            if in_gap:
                x_ranges.append((start_gap, len(c_masked_seq) - 1))

            chain_design_indices = []
            for start_idx, end_idx in x_ranges:
                r_start = start_idx + 1
                r_end = end_idx + 1
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
            # Calculate relative path from output_yaml_dir to pdb_path
            rel_path = os.path.relpath(str(pdb_path), str(output_yaml_dir))

            ab_entity = {
                "file": {
                    "path": rel_path,  # Relative path from YAML to PDB file
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

        # --- Process Antigen(s) ---
        if antigen_chains:
            for ag_chain in antigen_chains:
                ag_mapping = mapping

                # Compute epitopes
                epitopes = get_epitope_residues(pdb_path, pdb_path, h_chain or 'H', ag_chain)

                ag_include = [{"chain": {"id": ag_mapping.get(ag_chain, ag_chain)}}]
                ag_groups = [{"group": {"id": ag_mapping.get(ag_chain, ag_chain), "visibility": 2}}]
                ag_binding = []

                if epitopes:
                    range_str = format_residue_ranges(epitopes)
                    if range_str:
                        ag_binding.append({"chain": {"id": ag_mapping.get(ag_chain, ag_chain), "binding": range_str}})

                # Use same relative path as antibody
                ag_entity = {
                    "file": {
                        "path": rel_path,  # Same PDB file, same relative path
                        "include": ag_include,
                        "structure_groups": ag_groups,
                    }
                }
                if ag_binding:
                    ag_entity["file"]["binding_types"] = ag_binding

                yaml_entities.append(ag_entity)

        # --- Write YAML ---
        final_yaml = {"entities": yaml_entities}
        out_yaml_path = output_yaml_dir / f"{entry}.yaml"
        with open(out_yaml_path, 'w') as f:
            yaml.dump(final_yaml, f, sort_keys=False)

        logger.info(f"Generated YAML for {entry} -> {out_yaml_path}")
        processed_count += 1

    logger.info(f"\nSummary:")
    logger.info(f"  Processed: {processed_count}")
    logger.info(f"  Skipped: {skipped_count}")


if __name__ == "__main__":
    main()
