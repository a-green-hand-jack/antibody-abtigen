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
        output_dir: Optional[Path] = None,
        sabdab_summary_path: Optional[Path] = None
    ):
        """
        Initialize cleaner.

        Args:
            remove_water: Remove water molecules (HOH)
            remove_hetatm: Remove HETATM records (ligands, non-standard residues)
            output_dir: Directory to save cleaned CIF files (if None, use temp)
            sabdab_summary_path: Path to SAbDab summary TSV for chain type info
        """
        self.remove_water = remove_water
        self.remove_hetatm = remove_hetatm
        self.output_dir = output_dir
        self.sabdab_summary_path = sabdab_summary_path
        self._sabdab_cache = None

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

            # Step 2: Classify chains as antigen/antibody (using SAbDab metadata if available)
            chain_types = self._classify_chains(chain_info, pdb_id)

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

    def _load_sabdab_cache(self):
        """Load SAbDab summary file into cache."""
        if self._sabdab_cache is not None:
            return

        if self.sabdab_summary_path is None or not self.sabdab_summary_path.exists():
            self._sabdab_cache = {}
            return

        import csv
        self._sabdab_cache = {}

        with open(self.sabdab_summary_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pdb_id = row['pdb'].upper()
                self._sabdab_cache[pdb_id] = {
                    'Hchain': row.get('Hchain', ''),
                    'Lchain': row.get('Lchain', ''),
                    'antigen_chain': row.get('antigen_chain', ''),
                }

    def _classify_chains(
        self, chain_info: Dict[str, Dict], pdb_id: str
    ) -> Dict[str, str]:
        """
        Classify chains as antigen, antibody heavy, or antibody light.

        Uses SAbDab metadata if available, otherwise falls back to heuristics.

        Args:
            chain_info: Chain information dict
            pdb_id: PDB ID for SAbDab lookup

        Returns:
            Dict mapping chain_id → 'antigen' | 'antibody_heavy' | 'antibody_light'
        """
        chain_types = {}

        # Try SAbDab metadata first
        self._load_sabdab_cache()
        sabdab_info = self._sabdab_cache.get(pdb_id.upper(), {})

        hchain = sabdab_info.get('Hchain', '')
        lchain = sabdab_info.get('Lchain', '')
        antigen_chains = sabdab_info.get('antigen_chain', '').split('|')

        # Mark known chains from SAbDab
        for chain_id in chain_info.keys():
            if chain_id == hchain:
                chain_types[chain_id] = 'antibody_heavy'
            elif chain_id == lchain:
                chain_types[chain_id] = 'antibody_light'
            elif chain_id in antigen_chains:
                chain_types[chain_id] = 'antigen'

        # For unknown chains, use heuristics as fallback
        for chain_id, info in chain_info.items():
            if chain_id in chain_types:
                continue  # Already classified

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
        Create chain mappings (NO renaming - uses original PDB chain IDs).

        For epitope-centric pipeline, we don't rename chains because:
        1. Preserves original PDB structure integrity
        2. Avoids PyMOL visualization issues
        3. Easier to trace back to original PDB files

        The standardized_chain_id field is set to the original_chain_id.
        Chain type information is preserved in the chain_type field.

        Also builds auth_seq_id → 0-based index mapping.

        Args:
            pdb_id: PDB identifier
            chain_info: Chain information
            chain_types: Chain type classifications

        Returns:
            List of ChainMapping objects
        """
        mappings = []

        for chain_id in sorted(chain_info.keys()):
            chain_type = chain_types[chain_id]

            # Use original chain ID as standardized ID (no renaming)
            mappings.append(self._create_mapping(
                pdb_id, chain_id, chain_id, chain_type, chain_info
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
        Clean CIF file by filtering chains and removing water/HETATM.

        IMPORTANT: Does NOT rename chain IDs - uses original PDB chain IDs.
        Chain ID standardization is tracked in ChainMapping but not applied to files.

        Uses existing clean_structure() from structure.py which is battle-tested.

        Args:
            input_cif: Input CIF path
            pdb_id: PDB ID
            chain_mappings: Chain mappings (original_chain_id used)

        Returns:
            Path to cleaned CIF file
        """
        from ..structure import clean_structure as clean_cif

        # Determine output path
        if self.output_dir:
            output_path = self.output_dir / f"{pdb_id}_cleaned.cif"
            output_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            output_path = input_cif.parent / f"{pdb_id}_cleaned.cif"

        # Get original chain IDs to keep (no renaming)
        keep_chains = [m.original_chain_id for m in chain_mappings]

        # Use existing clean_structure function (battle-tested)
        clean_cif(
            input_cif=str(input_cif),
            output_cif=str(output_path),
            keep_chains=keep_chains,
            remove_water=self.remove_water,
            remove_hetatm=self.remove_hetatm
        )

        return output_path
