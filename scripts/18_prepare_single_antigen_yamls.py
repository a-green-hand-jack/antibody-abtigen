import argparse
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

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
    parser = argparse.ArgumentParser(description="Generate single-antigen YAML files for BoltzGen.")
    parser.add_argument("--input_dir", type=str, default="data/cleaned_split", help="Input directory containing paired antibody-antigen CIF files.")
    parser.add_argument("--output_cif_dir", type=str, default="data/antibody_antigen_single", help="Output directory for CIF files.")
    parser.add_argument("--output_yaml_dir", type=str, default="data/yamls/boltzgen_like/single_antigen", help="Output directory for generated YAML files.")
    parser.add_argument("--cdr_csv", type=str, default="data/decuplication/summary-decuplication-distance_threshold_9.csv", help="Path to CDR summary CSV.")
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
    """Maps Auth chain IDs to Label chain IDs from CIF file."""
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
    ns = gemmi.NeighborSearch(ab_struc[0], ab_struc.cell, dist_threshold + 2).populate()

    epitope_residues = {}

    # Iterate over Antigen atoms
    for model in ag_struc:
        for chain in model:
            chain_id = chain.name  # Auth ID
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


def get_structure_sequence(cif_path: Path, chain_id: str) -> str:
    """
    Extracts the sequence from the structure.
    Returns: sequence_string
    """
    struc = load_structure(cif_path)
    if not struc:
        return ""

    seq_str = ""

    for model in struc:
        for chain in model:
            if chain.name == chain_id:
                for res in chain:
                    try:
                        info = gemmi.find_tabulated_residue(res.name)
                        if info.is_amino_acid():
                            code = info.one_letter_code
                            seq_str += code
                    except:
                        pass
        break

    return seq_str


def main():
    args = parse_arguments()

    input_dir = Path(args.input_dir)
    output_cif_dir = Path(args.output_cif_dir)
    output_yaml_dir = Path(args.output_yaml_dir)

    setup_directories([output_cif_dir, output_yaml_dir])

    # Load CDR information
    logger.info(f"Loading CDR information from {args.cdr_csv}")
    cdr_df = pd.read_csv(args.cdr_csv)

    # Find all sample directories (containing paired CIFs)
    # Structure: data/cleaned_split/[cluster_id]/[pdb_chain]/
    samples = []
    for cluster_dir in input_dir.iterdir():
        if cluster_dir.is_dir():
            for sample_dir in cluster_dir.iterdir():
                if sample_dir.is_dir():
                    samples.append(sample_dir)

    samples = sorted(samples)
    logger.info(f"Found {len(samples)} sample directories.")

    processed_count = 0
    skipped_count = 0

    for sample_dir in samples:
        # Each sample_dir is like: data/cleaned_split/1793/7u0p_C/
        # Extract PDB code from directory name (e.g., 7u0p from 7u0p_C)
        pdb_chain = sample_dir.name
        pdb_code = pdb_chain.split('_')[0].lower()

        ab_cif = sample_dir / "antibody.cif"
        ag_cif = sample_dir / "antigen.cif"

        # Check if both files exist
        if not ab_cif.exists() or not ag_cif.exists():
            logger.warning(f"Skipping {pdb_chain}: Missing antibody.cif or antigen.cif")
            skipped_count += 1
            continue

        # Load structures
        ab_struc = load_structure(ab_cif)
        ag_struc = load_structure(ag_cif)

        if not ab_struc or not ag_struc:
            logger.warning(f"Skipping {pdb_chain}: Failed to load structures")
            skipped_count += 1
            continue

        # Get CDR information for this PDB
        cdr_rows = cdr_df[cdr_df['pdb'] == pdb_code]
        if cdr_rows.empty:
            logger.warning(f"Skipping {pdb_chain}: No CDR information found for {pdb_code}")
            skipped_count += 1
            continue

        cdr_info = cdr_rows.iloc[0]

        # --- Process Antibody ---
        ab_mapping = get_auth_to_label_mapping(ab_cif)

        yaml_entities = []

        # Antibody entity
        ab_include = []
        ab_design = []
        ab_structure_groups = []
        ab_exclude = []
        ab_design_insertions = []
        ab_reset = []

        chains_to_process = []
        if 'H_chain_id' in cdr_info and pd.notna(cdr_info['H_chain_id']):
            chains_to_process.append(('H', cdr_info['H_chain_id'], cdr_info.get('H_chain_seq'), cdr_info.get('H_chain_masked_seq')))
        if 'L_chain_id' in cdr_info and pd.notna(cdr_info['L_chain_id']):
            chains_to_process.append(('L', cdr_info['L_chain_id'], cdr_info.get('L_chain_seq'), cdr_info.get('L_chain_masked_seq')))

        for c_type, c_id_auth, c_seq, c_masked_seq in chains_to_process:
            if not c_id_auth or not c_masked_seq or not isinstance(c_masked_seq, str):
                continue

            # Convert Auth ID to Label ID
            c_id_label = ab_mapping.get(c_id_auth, c_id_auth)
            ab_include.append({"chain": {"id": c_id_label}})

            # Get masked residue ranges from CSV
            # Parse the masked sequence to find 'X' regions
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

            # For single antigen, we use CSV indices directly (they should match the structure)
            chain_design_indices = []
            for start_idx, end_idx in x_ranges:
                # Convert 0-based CSV indices to 1-based residue numbers
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
            # Save antibody CIF
            dst_ab_name = f"{pdb_chain}_antibody.cif"
            dst_ab_path = output_cif_dir / dst_ab_name
            ab_struc.make_mmcif_document().write_file(str(dst_ab_path))

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

        # --- Process Antigen ---
        ag_mapping = get_auth_to_label_mapping(ag_cif)

        # Compute epitopes
        epitopes = get_epitope_residues(ab_cif, ag_cif)

        # Get antigen chains
        ag_chains_auth = [c.name for c in ag_struc[0]]
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

        # Save antigen CIF
        dst_ag_name = f"{pdb_chain}_antigen.cif"
        dst_ag_path = output_cif_dir / dst_ag_name
        ag_struc.make_mmcif_document().write_file(str(dst_ag_path))

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

        # --- Write YAML ---
        final_yaml = {"entities": yaml_entities}
        out_yaml_path = output_yaml_dir / f"{pdb_chain}.yaml"
        with open(out_yaml_path, 'w') as f:
            yaml.dump(final_yaml, f, sort_keys=False)

        logger.info(f"Generated YAML for {pdb_chain} -> {out_yaml_path}")
        processed_count += 1

    logger.info(f"Processed {processed_count} samples, skipped {skipped_count}")


if __name__ == "__main__":
    main()
