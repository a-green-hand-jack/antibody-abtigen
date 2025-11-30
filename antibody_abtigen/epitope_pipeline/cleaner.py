"""
Structure cleaning module for the epitope-centric pipeline.

This module implements Phase 1: Structure Cleaning and Standardization.

Key tasks:
1. Parse mmCIF files using Gemmi
2. Filter for antigen + antibody chains
3. Standardize chain IDs (A,B,... for antigens; H,L for antibody)
4. Build auth_seq_id mapping for ESM-2 index conversion
5. Remove water and HETATM records

Reuses logic from:
- antibody_abtigen.structure.clean_structure()
- Adds chain ID standardization and index mapping
"""

import os
from pathlib import Path
from typing import Iterator, List, Dict, Optional
from collections import defaultdict

import gemmi
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

from .core import (
    StructureCleaner,
    CleanedStructure,
    ChainMapping,
    StructureCleaningError,
    IndexMappingError,
)


class GemmiStructureCleaner(StructureCleaner):
    """
    Structure cleaner using Gemmi library.

    This implementation:
    - Parses mmCIF files using Gemmi (fast and spec-compliant)
    - Identifies antigen and antibody chains
    - Standardizes chain IDs (antigen: A,B,C...; antibody: H,L)
    - Builds auth_seq_id → 0-based index mapping
    - Removes water (HOH) and HETATM records
    - Preserves mmCIF metadata categories

    Usage:
        >>> cleaner = GemmiStructureCleaner()
        >>> result = cleaner.clean_structure(Path("7K8T.cif"))
        >>> print(result.pdb_id)
        '7K8T'
        >>> print([m.standardized_chain_id for m in result.chain_mappings])
        ['A', 'B', 'H', 'L']
    """

    def __init__(
        self,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        output_dir: Optional[Path] = None
    ):
        """
        Initialize cleaner.

        Args:
            remove_water: Remove water molecules (HOH)
            remove_hetatm: Remove HETATM records (ligands, non-standard residues)
            output_dir: Directory to save cleaned CIF files (if None, use temp)
        """
        self.remove_water = remove_water
        self.remove_hetatm = remove_hetatm
        self.output_dir = output_dir

    def clean_structure(self, cif_path: Path) -> CleanedStructure:
        """
        Clean a single CIF structure.

        Args:
            cif_path: Path to input mmCIF file

        Returns:
            CleanedStructure with standardized chain IDs and mappings

        Raises:
            StructureCleaningError: If cleaning fails
        """
        if not cif_path.exists():
            raise StructureCleaningError(f"CIF file not found: {cif_path}")

        try:
            pdb_id = cif_path.stem.upper()

            # Step 1: Parse CIF with Biopython to get sequences
            chain_info = self._extract_chain_info(cif_path)

            # Step 2: Classify chains as antigen/antibody
            chain_types = self._classify_chains(chain_info)

            # Step 3: Create standardized chain mappings
            chain_mappings = self._create_chain_mappings(
                pdb_id, chain_info, chain_types
            )

            # Step 4: Clean CIF with Gemmi and save
            output_path = self._clean_and_save_cif(
                cif_path, pdb_id, chain_mappings
            )

            return CleanedStructure(
                pdb_id=pdb_id,
                file_path=output_path,
                chain_mappings=chain_mappings,
                antigen_cluster_id=None  # Assigned later in pipeline
            )

        except Exception as e:
            raise StructureCleaningError(
                f"Failed to clean {cif_path}: {str(e)}"
            ) from e

    def batch_clean(
        self, cif_paths: List[Path]
    ) -> Iterator[CleanedStructure]:
        """
        Clean multiple structures (generator for memory efficiency).

        Args:
            cif_paths: List of input CIF paths

        Yields:
            CleanedStructure objects
        """
        for cif_path in cif_paths:
            try:
                yield self.clean_structure(cif_path)
            except StructureCleaningError as e:
                # Log error but continue processing
                print(f"Warning: {e}")
                continue

    def _extract_chain_info(self, cif_path: Path) -> Dict[str, Dict]:
        """
        Extract chain sequences and metadata from CIF.

        Args:
            cif_path: Path to CIF file

        Returns:
            Dict mapping chain_id → {sequence, length, residue_nums}
        """
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("tmp", str(cif_path))
        model = structure[0]

        chain_info = {}
        for chain in model:
            chain_id = chain.get_id()

            # Extract residues (standard amino acids only)
            residues = [r for r in chain if r.get_id()[0] == ' ']

            # Get sequence
            sequence = ""
            auth_seq_ids = []
            for res in residues:
                try:
                    aa = seq1(res.get_resname())
                    if aa != 'X':  # Skip unknown residues
                        sequence += aa
                        auth_seq_ids.append(res.get_id()[1])  # auth_seq_id
                except (KeyError, ValueError):
                    pass

            if sequence:  # Only include chains with valid sequence
                chain_info[chain_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'auth_seq_ids': auth_seq_ids
                }

        return chain_info

    def _classify_chains(
        self, chain_info: Dict[str, Dict]
    ) -> Dict[str, str]:
        """
        Classify chains as antigen, antibody heavy, or antibody light.

        Uses SAbDab heuristics:
        - Heavy chain: 250-500 residues (typical VH-CH1-CH2-CH3)
        - Light chain: 200-250 residues (typical VL-CL)
        - Antigen: Everything else

        Args:
            chain_info: Chain information dict

        Returns:
            Dict mapping chain_id → 'antigen' | 'antibody_heavy' | 'antibody_light'
        """
        chain_types = {}

        for chain_id, info in chain_info.items():
            length = info['length']

            # Heavy chain: typically 250-500 residues
            if 250 <= length <= 500:
                chain_types[chain_id] = 'antibody_heavy'
            # Light chain: typically 200-250 residues
            elif 200 <= length <= 250:
                chain_types[chain_id] = 'antibody_light'
            # Everything else is antigen
            else:
                chain_types[chain_id] = 'antigen'

        return chain_types

    def _create_chain_mappings(
        self,
        pdb_id: str,
        chain_info: Dict[str, Dict],
        chain_types: Dict[str, str]
    ) -> List[ChainMapping]:
        """
        Create standardized chain mappings.

        Standardization rules:
        - Antigen chains: A, B, C, D, ... (alphabetical order by original ID)
        - Antibody heavy: H
        - Antibody light: L

        Also builds auth_seq_id → 0-based index mapping.

        Args:
            pdb_id: PDB identifier
            chain_info: Chain information
            chain_types: Chain type classifications

        Returns:
            List of ChainMapping objects
        """
        # Separate chains by type
        antigen_chains = []
        heavy_chains = []
        light_chains = []

        for chain_id in sorted(chain_info.keys()):
            chain_type = chain_types[chain_id]
            if chain_type == 'antigen':
                antigen_chains.append(chain_id)
            elif chain_type == 'antibody_heavy':
                heavy_chains.append(chain_id)
            elif chain_type == 'antibody_light':
                light_chains.append(chain_id)

        # Create mappings
        mappings = []

        # Antigen chains: A, B, C, ...
        for idx, original_id in enumerate(antigen_chains):
            standardized_id = chr(ord('A') + idx)  # A, B, C, ...
            mappings.append(self._create_mapping(
                pdb_id, original_id, standardized_id, 'antigen', chain_info
            ))

        # Antibody heavy: H (take first heavy chain if multiple)
        if heavy_chains:
            mappings.append(self._create_mapping(
                pdb_id, heavy_chains[0], 'H', 'antibody_heavy', chain_info
            ))

        # Antibody light: L (take first light chain if multiple)
        if light_chains:
            mappings.append(self._create_mapping(
                pdb_id, light_chains[0], 'L', 'antibody_light', chain_info
            ))

        return mappings

    def _create_mapping(
        self,
        pdb_id: str,
        original_id: str,
        standardized_id: str,
        chain_type: str,
        chain_info: Dict[str, Dict]
    ) -> ChainMapping:
        """
        Create a single ChainMapping with index map.

        Args:
            pdb_id: PDB ID
            original_id: Original chain ID from CIF
            standardized_id: Standardized chain ID (A, B, H, L, etc.)
            chain_type: 'antigen' | 'antibody_heavy' | 'antibody_light'
            chain_info: Chain information dict

        Returns:
            ChainMapping with auth_seq_id_map populated
        """
        info = chain_info[original_id]
        sequence = info['sequence']
        auth_seq_ids = info['auth_seq_ids']

        # Build mapping: 0-based index → auth_seq_id
        auth_seq_id_map = {
            idx: auth_id
            for idx, auth_id in enumerate(auth_seq_ids)
        }

        return ChainMapping(
            pdb_id=pdb_id,
            original_chain_id=original_id,
            standardized_chain_id=standardized_id,
            chain_type=chain_type,
            sequence=sequence,
            auth_seq_id_map=auth_seq_id_map
        )

    def _clean_and_save_cif(
        self,
        input_cif: Path,
        pdb_id: str,
        chain_mappings: List[ChainMapping]
    ) -> Path:
        """
        Clean CIF file and save with standardized chain IDs.

        Uses Gemmi to:
        - Filter chains
        - Remove water/HETATM
        - Rename chains to standardized IDs
        - Preserve metadata

        Args:
            input_cif: Input CIF path
            pdb_id: PDB ID
            chain_mappings: Chain mappings with standardization

        Returns:
            Path to cleaned CIF file
        """
        # Determine output path
        if self.output_dir:
            output_path = self.output_dir / f"{pdb_id}_cleaned.cif"
            output_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            output_path = input_cif.parent / f"{pdb_id}_cleaned.cif"

        # Build chain ID mapping
        original_to_std = {
            m.original_chain_id: m.standardized_chain_id
            for m in chain_mappings
        }
        keep_chains = set(original_to_std.keys())

        # Load CIF with Gemmi
        doc = gemmi.cif.read_file(str(input_cif))
        block = doc.sole_block()

        # Find atom_site table
        table = block.find_mmcif_category("_atom_site")
        if not table:
            # No atoms, just write original
            doc.write_file(str(output_path))
            return output_path

        tags = list(table.tags)

        def get_idx(tag: str) -> int:
            """Get column index for a tag."""
            full_tag = tag if tag.startswith("_") else f"_atom_site.{tag}"
            try:
                return tags.index(full_tag)
            except ValueError:
                return -1

        # Get column indices
        idx_auth = get_idx("_atom_site.auth_asym_id")
        idx_label = get_idx("_atom_site.label_asym_id")
        idx_comp = get_idx("_atom_site.label_comp_id")
        idx_group = get_idx("_atom_site.group_PDB")

        # Filter and rename rows
        rows_to_remove = []
        for i, row in enumerate(table):
            auth_chain = row[idx_auth] if idx_auth >= 0 else ""
            label_chain = row[idx_label] if idx_label >= 0 else ""
            comp_id = row[idx_comp] if idx_comp >= 0 else ""
            group = row[idx_group] if idx_group >= 0 else ""

            keep = True

            # Filter by chain
            if auth_chain not in keep_chains:
                keep = False

            # Remove water
            if keep and self.remove_water and comp_id == "HOH":
                keep = False

            # Remove HETATM
            if keep and self.remove_hetatm and group == "HETATM":
                keep = False

            if not keep:
                rows_to_remove.append(i)
            else:
                # Rename chain to standardized ID
                if auth_chain in original_to_std:
                    new_chain = original_to_std[auth_chain]
                    if idx_auth >= 0:
                        row[idx_auth] = new_chain
                    if idx_label >= 0:
                        row[idx_label] = new_chain

        # Remove filtered rows (reverse order to preserve indices)
        for offset, row_idx in enumerate(rows_to_remove):
            table.remove_row(row_idx - offset)

        # Check for empty loops
        block.check_empty_loops("_atom_site")

        # Write output
        doc.write_file(str(output_path))

        return output_path
