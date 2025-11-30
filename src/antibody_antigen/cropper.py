"""
AntiGenCropper: Extracts the local environment (pocket) around an epitope.
"""

from typing import List, Dict, Set, Optional, Tuple
import numpy as np
import gemmi

class AntiGenCropper:
    """
    Implements the Affinity Cropper algorithm (Algorithm 3) to extract
    antigen pockets around defined epitope residues.
    """
    def __init__(
        self, 
        neighborhood_size: int = 10, 
        max_tokens: int = 256, 
        max_protein_tokens: int = 200
    ):
        """
        Args:
            neighborhood_size: Number of residues to include in the window around each epitope residue.
            max_tokens: Maximum total tokens (atoms or residues) - reserve for future use.
            max_protein_tokens: Maximum number of protein residues to keep in the crop.
        """
        self.neighborhood_size = neighborhood_size
        self.max_tokens = max_tokens
        self.max_protein_tokens = max_protein_tokens

    def crop(
        self, 
        structure: gemmi.Structure, 
        epitope_residues_map: Dict[str, List[gemmi.Residue]]
    ) -> Dict[str, np.ndarray]:
        """
        Extracts the pocket environment around the provided epitope residues.

        Args:
            structure: The full gemmi Structure (or Model).
            epitope_residues_map: Dict mapping chain ID to list of epitope residues in that chain.

        Returns:
            A dictionary containing numpy arrays of the cropped structure data.
        """
        if not epitope_residues_map:
            return {}

        # 1. Use the first model if structure is provided
        if isinstance(structure, gemmi.Structure):
            model = structure[0]
        else:
            model = structure
        
        # epitope_residues_map is already organized by chain: {chain_name: [res...]}
        epitope_by_chain = epitope_residues_map
            
        cropped_residues_set: Set[str] = set() # Stores unique ID "chain_seqid"
        cropped_residues_list: List[Tuple[str, gemmi.Residue]] = []
        
        # 2. Process each chain
        for chain_name, epi_res_list in epitope_by_chain.items():
            try:
                chain = model[chain_name]
            except LookupError:
                continue
                
            # Get all residues in chain to allow index-based window expansion
            all_chain_res = list(chain)
            # Map seqid to index in all_chain_res
            # We use str(seqid) to handle insertion codes if necessary
            seqid_to_idx = {str(r.seqid): i for i, r in enumerate(all_chain_res)}
            
            # Identify indices of epitope residues
            epi_indices = []
            for r in epi_res_list:
                sid = str(r.seqid)
                if sid in seqid_to_idx:
                    epi_indices.append(seqid_to_idx[sid])
            
            epi_indices.sort()
            
            # Algorithm 3 Logic: Expand window for each epitope residue
            # "res_idx_sorted" -> in our case, we iterate through the provided epitope residues
            # (which are usually sorted by distance in the caller, but here we iterate all)
            
            chain_cropped_indices = set()
            
            for core_idx in epi_indices:
                min_idx = core_idx
                max_idx = core_idx
                current_window = {core_idx}
                
                # Expand
                # Note: The pseudo-code expands until `len(res_tokens) < neighborhood_size`
                # `res_tokens` seems to be the tokens FOR THE CURRENT RESIDUE's window.
                while len(current_window) < self.neighborhood_size:
                    expanded = False
                    # Try expand left
                    if min_idx > 0:
                        min_idx -= 1
                        current_window.add(min_idx)
                        expanded = True
                    
                    # Try expand right (only if we still need more and can expand)
                    if len(current_window) < self.neighborhood_size and max_idx < len(all_chain_res) - 1:
                        max_idx += 1
                        current_window.add(max_idx)
                        expanded = True
                        
                    if not expanded:
                        break
                
                chain_cropped_indices.update(current_window)
            
            # Add to global list
            sorted_indices = sorted(list(chain_cropped_indices))
            for idx in sorted_indices:
                res = all_chain_res[idx]
                unique_id = f"{chain_name}_{res.seqid}"
                if unique_id not in cropped_residues_set:
                    # Check Global Token Limit (Residue Level)
                    if len(cropped_residues_list) >= self.max_protein_tokens:
                        break 
                    
                    cropped_residues_set.add(unique_id)
                    cropped_residues_list.append((chain_name, res))
            
            if len(cropped_residues_list) >= self.max_protein_tokens:
                break

        # 3. Extract Atomic Data
        # We want to return data structured for easy ML usage or storage
        coords = []
        atom_types = []
        residue_indices = [] # 0..N_res-1
        chain_ids = []
        res_seq_nums = []
        res_names = []
        
        # Sort residues? 
        # It's better to keep them grouped by chain, then by seq id.
        # They are currently appended in chain loop order, and sorted by index within chain.
        # That is a good natural order.
        
        for local_res_idx, (c_id, res) in enumerate(cropped_residues_list):
            r_seq = res.seqid.num
            r_name = res.name
            
            for atom in res:
                # Skip hydrogens? Usually yes for rough geometry, but gemmi keeps them if present.
                # Let's keep all atoms for now, user didn't specify to strip H.
                # Standard PDB processing often strips H. 
                # If `element` is H or D, maybe skip? 
                # Let's keep them, it's safer to filter later than miss them now.
                
                coords.append(atom.pos.tolist())
                atom_types.append(atom.element.name)
                residue_indices.append(local_res_idx)
                chain_ids.append(c_id)
                res_seq_nums.append(r_seq)
                res_names.append(r_name)
                
        if not coords:
            return {}

        return {
            "coords": np.array(coords, dtype=np.float32),
            "atom_types": np.array(atom_types, dtype=str),
            "residue_indices": np.array(residue_indices, dtype=np.int32),
            "chain_ids": np.array(chain_ids, dtype=str),
            "res_seq_nums": np.array(res_seq_nums, dtype=np.int32),
            "res_names": np.array(res_names, dtype=str)
        }
