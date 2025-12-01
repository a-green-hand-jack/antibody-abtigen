import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import collections

import gemmi
import numpy as np
import pandas as pd
import yaml

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants for paths
CLEANED_SPLIT_DIR = Path("data/cleaned_split")
CLUSTER_SUMMARY_CSV = Path("data/cluster/epitope_cluster_summary.csv")
SABDAB_SUMMARY_TSV = Path("data/raw_data/meta/sabdab_summary_all.tsv")
POCKETS_NPZ_DIR = Path("data/extract_epitope/pockets")
OUTPUT_YAML_DIR = Path("data/yamls/boltzgen_like/multi_antigen")

def parse_cif_chains(cif_path: Path) -> List[str]:
    """
    Parses a CIF file to extract unique chain IDs.
    Uses gemmi.cif.Block.find_mmcif_category to robustly find categories.
    Prioritizes:
    1. _entity_poly.pdbx_strand_id
    2. _atom_site.label_asym_id
    3. _atom_site.auth_asym_id
    """
    try:
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        
        chain_ids = set()

        # 1. Try _entity_poly.pdbx_strand_id
        # This is the canonical place for chain IDs in mmCIF for polymers
        poly_cat = block.find_mmcif_category('_entity_poly.')
        if poly_cat:
            try:
                for val in poly_cat.find_column('pdbx_strand_id'):
                    # Chains can be comma-separated in this field
                    for chain_id_part in val.replace("'", "").replace('"', '').split(','):
                        chain_id = chain_id_part.strip()
                        if chain_id and chain_id not in ['?', '.']:
                            chain_ids.add(chain_id)
            except ValueError:
                pass # Column not found

        if chain_ids:
            return sorted(list(chain_ids))

        # 2. Try _atom_site.label_asym_id (canonical for atoms) or auth_asym_id (author provided)
        atom_cat = block.find_mmcif_category('_atom_site.')
        if atom_cat:
            # Try label_asym_id
            try:
                for val in atom_cat.find_column('label_asym_id'):
                    chain_id = val.strip()
                    if chain_id and chain_id not in ['?', '.']:
                        chain_ids.add(chain_id)
            except ValueError:
                pass
            
            if chain_ids:
                return sorted(list(chain_ids))

            # Try auth_asym_id
            try:
                for val in atom_cat.find_column('auth_asym_id'):
                    chain_id = val.strip()
                    if chain_id and chain_id not in ['?', '.']:
                        chain_ids.add(chain_id)
            except ValueError:
                pass

        return sorted(list(chain_ids))

    except Exception as e:
        logger.warning(f"Could not parse chains from {cif_path}: {e}")
        return []

def load_epitope_residues(npz_path: Path) -> Dict[str, List[int]]:
    """
    Loads residue sequence numbers for epitope from an NPZ file.
    """
    if not npz_path.exists():
        logger.warning(f"NPZ file not found: {npz_path}")
        return {}
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        chain_ids = data['chain_ids']
        res_seq_nums = data['res_seq_nums']

        epitope_map: Dict[str, List[int]] = collections.defaultdict(list)
        for chain_id, res_seq_num in zip(chain_ids, res_seq_nums):
            epitope_map[str(chain_id)].append(int(res_seq_num))
        
        # Sort and unique residue numbers for each chain
        for chain_id in epitope_map:
            epitope_map[chain_id] = sorted(list(set(epitope_map[chain_id])))

        return dict(epitope_map)
    except Exception as e:
        logger.error(f"Failed to load epitope residues from {npz_path}: {e}")
        return {}

def get_reference_antibody(
    cluster_members_df: pd.DataFrame, 
    sabdab_summary_df: pd.DataFrame
) -> Optional[Tuple[str, pd.Series]]:
    """
    Selects a reference antibody epitope_id based on species preference (Human > Mouse > First).
    Returns (epitope_id, corresponding row from cluster_members_df).
    """
    
    potential_refs = cluster_members_df.copy()
    # Ensure PDB IDs are lowercase for matching sabdab_summary
    potential_refs['pdb_id_lower'] = potential_refs['pdb_id'].str.lower()

    # Merge with sabdab summary to get species info
    # We need H_chain, L_chain for correct matching as PDB ID can have multiple entries
    merged_df = pd.merge(
        potential_refs,
        sabdab_summary_df,
        left_on=['pdb_id_lower', 'H_chain', 'L_chain'],
        right_on=['pdb', 'Hchain', 'Lchain'],
        how='left',
        suffixes=('_cluster', '_sabdab')
    )
    
    # Filter out entries where species info isn't available
    merged_df = merged_df.dropna(subset=['heavy_species', 'light_species'])

    if merged_df.empty:
        logger.warning(f"No species information found for cluster members.")
        # Fallback to just taking the first one if no species info is found
        if not cluster_members_df.empty:
            first_row = cluster_members_df.iloc[0]
            return first_row['epitope_id'], first_row
        return None, None

    # Prioritize Human
    human_abs = merged_df[
        (merged_df['heavy_species'].str.contains('Homo sapiens', case=False, na=False)) |
        (merged_df['light_species'].str.contains('Homo sapiens', case=False, na=False))
    ]
    if not human_abs.empty:
        # Choose the first one that is a representative if available
        human_representative = human_abs[human_abs['is_representative'] == True]
        if not human_representative.empty:
            ref_row = human_representative.iloc[0]
        else:
            ref_row = human_abs.iloc[0]
        return ref_row['epitope_id'], ref_row

    # Then Mouse
    mouse_abs = merged_df[
        (merged_df['heavy_species'].str.contains('Mus musculus', case=False, na=False)) |
        (merged_df['light_species'].str.contains('Mus musculus', case=False, na=False))
    ]
    if not mouse_abs.empty:
        mouse_representative = mouse_abs[mouse_abs['is_representative'] == True]
        if not mouse_representative.empty:
            ref_row = mouse_representative.iloc[0]
        else:
            ref_row = mouse_abs.iloc[0]
        return ref_row['epitope_id'], ref_row
    
    # Fallback to the first available member if no specific species found
    if not cluster_members_df.empty:
        first_row = cluster_members_df.iloc[0]
        logger.warning(f"No Human or Mouse antibody found in cluster, using first entry: {first_row['epitope_id']}")
        return first_row['epitope_id'], first_row
        
    return None, None

def main():
    OUTPUT_YAML_DIR.mkdir(parents=True, exist_ok=True)

    # Load necessary metadata
    cluster_summary_df = pd.read_csv(CLUSTER_SUMMARY_CSV)
    sabdab_summary_df = pd.read_csv(SABDAB_SUMMARY_TSV, sep='\t', low_memory=False)
    
    # Ensure some sabdab columns are strings for merging
    sabdab_summary_df['pdb'] = sabdab_summary_df['pdb'].astype(str).str.lower()
    sabdab_summary_df['Hchain'] = sabdab_summary_df['Hchain'].fillna('').astype(str)
    sabdab_summary_df['Lchain'] = sabdab_summary_df['Lchain'].fillna('').astype(str)

    # Group cluster_summary by cluster_id
    for cluster_id_str in sorted(
        [d.name for d in CLEANED_SPLIT_DIR.iterdir() if d.is_dir()]
    ):
        try:
            cluster_id = int(cluster_id_str)
        except ValueError:
            logger.warning(f"Skipping non-integer cluster directory: {cluster_id_str}")
            continue

        cluster_dir_path = CLEANED_SPLIT_DIR / cluster_id_str
        
        logger.info(f"Processing cluster: {cluster_id}")

        # Get all members of the current cluster from the summary CSV
        cluster_members_df = cluster_summary_df[cluster_summary_df['cluster_id'] == cluster_id]
        
        if cluster_members_df.empty:
            logger.warning(f"No entries found in epitope_cluster_summary.csv for cluster_id {cluster_id}, skipping.")
            continue

        # Select the reference antibody
        ref_epitope_id, ref_row = get_reference_antibody(cluster_members_df, sabdab_summary_df)

        if ref_epitope_id is None:
            logger.error(f"Could not determine a reference antibody for cluster {cluster_id}, skipping.")
            continue
        
        logger.info(f"Reference antibody for cluster {cluster_id}: {ref_epitope_id}")
        
        # --- Construct Antibody Section ---
        antibody_yaml_data: Dict[str, Any] = {}
        ref_antibody_cif_path = cluster_dir_path / ref_epitope_id / "antibody.cif"
        
        if not ref_antibody_cif_path.exists():
            logger.error(f"Reference antibody CIF not found for {ref_epitope_id} at {ref_antibody_cif_path}, skipping cluster.")
            continue

        ab_chains_in_cif = parse_cif_chains(ref_antibody_cif_path)
        antibody_yaml_data["file"] = {
            "path": str(ref_antibody_cif_path),
            "include": [{"chain": {"id": c}} for c in ab_chains_in_cif]
        }
        antibody_yaml_data["chains"] = ab_chains_in_cif

        # --- Construct Antigens Section ---
        antigens_yaml_list: List[Dict[str, Any]] = []
        
        # Iterate through all member epitope_id directories within the cluster
        for member_epitope_dir in sorted(cluster_dir_path.iterdir()):
            if not member_epitope_dir.is_dir():
                continue
            
            member_epitope_id = member_epitope_dir.name
            
            member_antigen_cif_path = member_epitope_dir / "antigen.cif"
            if not member_antigen_cif_path.exists():
                logger.warning(f"Antigen CIF not found for {member_epitope_id} at {member_antigen_cif_path}, skipping this antigen.")
                continue

            ag_chains_in_cif = parse_cif_chains(member_antigen_cif_path)
            if not ag_chains_in_cif:
                logger.warning(f"No antigen chains found in {member_antigen_cif_path}, skipping this antigen.")
                continue
            
            npz_path = POCKETS_NPZ_DIR / f"{member_epitope_id}.npz"
            epitope_res_map = load_epitope_residues(npz_path)

            antigen_entry: Dict[str, Any] = {
                "name": member_epitope_id,
                "file": {
                    "path": str(member_antigen_cif_path),
                    "include": [{"chain": {"id": c}} for c in ag_chains_in_cif]
                }
            }
            
            epitope_list_for_yaml = []
            for chain, res_ids in epitope_res_map.items():
                # Only include epitope residues for chains actually present in the antigen CIF
                if chain in ag_chains_in_cif:
                    epitope_list_for_yaml.append({
                        "chain": chain,
                        "residue_ids": res_ids
                    })
            
            if epitope_list_for_yaml:
                antigen_entry["epitope"] = epitope_list_for_yaml
            else:
                logger.warning(f"No epitope residues found or matched for {member_epitope_id} in {member_antigen_cif_path}.")


            antigens_yaml_list.append(antigen_entry)
        
        if not antigens_yaml_list:
            logger.error(f"No valid antigens found for cluster {cluster_id}, skipping YAML generation.")
            continue

        # --- Assemble Full YAML ---
        full_yaml_content = {
            "multi_antigen": {
                "antibody": antibody_yaml_data,
                "antigens": antigens_yaml_list,
                "parallel_config": {
                    "aggregation_method": "kabsch_mean",
                    "antigen_weights": None,
                },
            }
            # "constraints" field is omitted as per discussion
        }

        output_yaml_path = OUTPUT_YAML_DIR / f"cluster_{cluster_id}.yaml"
        with open(output_yaml_path, 'w') as f:
            yaml.dump(full_yaml_content, f, sort_keys=False, indent=2)
        
        logger.info(f"Generated YAML for cluster {cluster_id} at {output_yaml_path}")

if __name__ == "__main__":
    main()
