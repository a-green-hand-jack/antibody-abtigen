import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import gemmi
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

# Add src to path to import internal modules
sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))

try:
    from antibody_antigen.cropper import AntiGenCropper
except ImportError:
    # Fallback if running from different context
    sys.path.append('src')
    from antibody_antigen.cropper import AntiGenCropper

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
DIST_THRESHOLD_EPITOPE = 4.5  # Angstroms

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract epitopes and affinity pockets from antibody-antigen complexes.")
    parser.add_argument(
        "--cif_dir",
        type=str,
        required=True,
        help="Directory containing input CIF files."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Path to the input CSV file (e.g., deduplicated summary)."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save output CSV and NPZ files."
    )
    parser.add_argument(
        "--min_residues",
        type=int,
        default=6,
        help="Minimum number of residues in an epitope to keep."
    )
    parser.add_argument(
        "--max_residues",
        type=int,
        default=60,
        help="Maximum number of residues in an epitope to keep."
    )
    parser.add_argument(
        "--neighborhood_size",
        type=int,
        default=10,
        help="Window size for affinity cropping (Algorithm 3)."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of parallel workers."
    )
    return parser.parse_args()

def get_structure_from_cif(cif_path: Path) -> Optional[gemmi.Structure]:
    """Reads a CIF file and returns the gemmi Structure (Model 1)."""
    try:
        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)
        return structure
    except Exception as e:
        # logger.debug(f"Failed to read CIF {cif_path}: {e}")
        return None

def extract_epitope_residues(
    structure: gemmi.Structure,
    h_chain_id: str,
    l_chain_id: str,
    ag_chain_ids: List[str],
    dist_threshold: float = DIST_THRESHOLD_EPITOPE
) -> Dict[str, List[gemmi.Residue]]:
    """
    Identifies antigen residues within distance threshold of antibody (H/L) chains.
    Returns a dictionary mapping antigen chain ID to list of epitope residues.
    """
    model = structure[0] # Assume first model
    
    # 1. Collect Antibody Chains for NeighborSearch
    # Create a temporary model containing only antibody chains to populate NeighborSearch
    ab_model = gemmi.Model("1")
    
    def add_chain_to_temp_model(chain_id):
        if not chain_id or chain_id == 'nan': return
        try:
            # Clone the chain to add to temp model
            original_chain = model[chain_id]
            ab_model.add_chain(original_chain)
        except LookupError:
            pass 

    add_chain_to_temp_model(h_chain_id)
    add_chain_to_temp_model(l_chain_id)
    
    if len(ab_model) == 0:
        return []

    # Build NeighborSearch for Antibody
    # Populating with the temporary antibody model ensures only antibody atoms are in the tree
    ns = gemmi.NeighborSearch(ab_model, structure.cell, max(5.0, dist_threshold)).populate(include_h=False)
    
    # 2. Find Antigen residues close to Antibody
    epitope_residues_map: Dict[str, List[gemmi.Residue]] = {}
    
    for ag_id in ag_chain_ids:
        try:
            ag_chain = model[ag_id]
        except LookupError:
            continue

        chain_epitope_list = []
        seen_res = set()

        for res in ag_chain:
            is_epitope = False
            for atom in res:
                # Check if any Antibody atom is within threshold
                marks = ns.find_atoms(atom.pos, radius=dist_threshold)
                if len(marks) > 0:
                    is_epitope = True
                    break # One atom is enough to mark residue
            
            if is_epitope:
                # Unique check using seqid within chain
                sid = str(res.seqid)
                if sid not in seen_res:
                    seen_res.add(sid)
                    chain_epitope_list.append(res)
        
        if chain_epitope_list:
            epitope_residues_map[ag_id] = chain_epitope_list

    return epitope_residues_map

def process_entry(
    row: pd.Series, 
    cif_dir: Path, 
    output_npz_dir: Path,
    cropper: AntiGenCropper,
    min_res: int,
    max_res: int
) -> Optional[Dict]:
    
    pdb_id = str(row['pdb'])
    
    # Attempt to find file
    cif_path = cif_dir / f"{pdb_id}.cif"
    if not cif_path.exists():
        cif_path = cif_dir / f"{pdb_id.lower()}.cif"
        if not cif_path.exists():
            return None

    structure = get_structure_from_cif(cif_path)
    if not structure:
        return None
        
    # Parse Chain IDs
    h_chain = str(row['H_chain_id']) if pd.notna(row['H_chain_id']) else ''
    l_chain = str(row['L_chain_id']) if pd.notna(row['L_chain_id']) else ''
    ag_chain_str = str(row['antigen_chain_id']) if pd.notna(row['antigen_chain_id']) else ''
    
    # Handle multi-chain antigen format "A|B" or "A"
    ag_chains = []
    if ag_chain_str and ag_chain_str != 'nan':
        # Clean separators and Python list artifacts
        clean_str = ag_chain_str.replace('|', ',').replace('[', '').replace(']', '').replace("'", "").replace('"', '')
        parts = [p.strip() for p in clean_str.split(',') if p.strip()]
        ag_chains = sorted(parts)
    
    if not ag_chains:
        return None

    # 1. Extract Epitope (Raw)
    # Treating all antigen chains as ONE system
    epitope_res_map = extract_epitope_residues(structure, h_chain, l_chain, ag_chains)
    
    # Calculate total residues
    num_res = sum(len(v) for v in epitope_res_map.values())
    if num_res < min_res or num_res > max_res:
        return None # Filtered out
            
    # 2. Affinity Cropping (Context)
    try:
        pocket_data = cropper.crop(structure, epitope_res_map)
    except Exception as e:
        logger.warning(f"Cropping failed for {pdb_id}: {e}")
        return None
    
    if not pocket_data:
        return None
        
    # 3. Save
    # ID Generation: PDB_AgChains
    # Normalize chain string for ID: "A_B"
    ag_chain_id_flat = "_".join(ag_chains)
    epitope_id = f"{pdb_id}_{ag_chain_id_flat}"
    
    npz_filename = f"{epitope_id}.npz"
    npz_path = output_npz_dir / npz_filename
    
    np.savez_compressed(npz_path, **pocket_data)
    
    # Return Metadata Record
    record = {
        'epitope_id': epitope_id,
        'pdb_id': pdb_id,
        'pdb_path': str(cif_path),
        'H_chain': h_chain,
        'L_chain': l_chain,
        'antigen_chains': ag_chain_str, # Keep original format
        'num_epitope_residues': num_res,
        'num_pocket_atoms': len(pocket_data['coords']),
        'pocket_npz_path': str(npz_path),
        'original_index': row.get('index_in_summary', -1)
    }
        
    return record

def main():
    args = parse_arguments()
    
    # Setup directories
    cif_dir = Path(args.cif_dir)
    output_dir = Path(args.output_dir)
    output_npz_dir = output_dir / "pockets"
    output_npz_dir.mkdir(parents=True, exist_ok=True)
    
    # Load Metadata
    if not Path(args.input_csv).exists():
        logger.error(f"Input CSV not found: {args.input_csv}")
        return

    df = pd.read_csv(args.input_csv)
    logger.info(f"Loaded {len(df)} entries from {args.input_csv}")
    
    # Initialize Cropper
    cropper = AntiGenCropper(neighborhood_size=args.neighborhood_size)
    
    # Parallel Processing
    # Note: Gemmi objects (Structure) are not easily picklable for multiprocessing 
    # if passed directly, but here we process per-row inside the worker function.
    # process_entry takes the row and paths, handles IO internally.
    
    results_list = Parallel(n_jobs=args.workers)(
        delayed(process_entry)(
            row, 
            cif_dir, 
            output_npz_dir, 
            cropper,
            args.min_residues,
            args.max_residues
        )
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing Epitopes")
    )
    
    # Filter Nones
    valid_results = [r for r in results_list if r is not None]
    
    # Save Summary
    if valid_results:
        result_df = pd.DataFrame(valid_results)
        output_csv_path = output_dir / "epitope_summary.csv"
        result_df.to_csv(output_csv_path, index=False)
        logger.info(f"Completed. Saved {len(result_df)} records to {output_csv_path}")
    else:
        logger.warning("No valid epitopes found matching criteria.")

if __name__ == "__main__":
    main()
