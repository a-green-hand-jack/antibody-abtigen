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
from typing import Iterator, List, Dict, Optional, Tuple
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum

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


class FilterReason(Enum):
    """Reasons for filtering out structures."""
    ACCEPTED = "accepted"
    NO_ANTIGEN = "no_antigen"
    PEPTIDE_ONLY = "peptide_only"
    HAPTEN_ONLY = "hapten_only"
    NO_PROTEIN_ANTIGEN = "no_protein_antigen"
    INSUFFICIENT_ANTIGEN_RESIDUES = "insufficient_antigen_residues"


@dataclass(frozen=True)
class FilterResult:
    """Result of structure filtering decision."""
    pdb_id: str
    accepted: bool
    reason: FilterReason
    antigen_chains: str  # e.g., "A | B"
    antigen_type: str    # e.g., "protein | peptide"
    num_antigen_chains: int
    total_antigen_residues: int
    num_antibody_chains: int
    details: str  # Additional human-readable details


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

    def should_process_structure(
        self,
        pdb_id: str,
        chain_info: Dict[str, Dict],
        chain_types: Dict[str, str]
    ) -> FilterResult:
        """
        Determine if structure should be processed based on antigen type.

        Filtering rules:
        1. ACCEPT: Has protein antigen (pure protein or protein+peptide)
        2. REJECT: No antigen (antigen_chain=NA)
        3. REJECT: Peptide-only antigen (no protein component)
        4. REJECT: Hapten-only antigen
        5. REJECT: Insufficient antigen residues (<10 residues total)

        Args:
            pdb_id: PDB identifier
            chain_info: Chain sequence information
            chain_types: Chain type classifications

        Returns:
            FilterResult with decision and details
        """
        # Load SAbDab metadata
        self._load_sabdab_cache()
        sabdab_info = self._sabdab_cache.get(pdb_id.upper(), {})

        antigen_chain_str = sabdab_info.get('antigen_chain', '')
        antigen_type_str = sabdab_info.get('antigen_type', '')

        # Count chains
        antigen_chains = [cid for cid, ctype in chain_types.items() if ctype == 'antigen']
        antibody_chains = [cid for cid, ctype in chain_types.items() if ctype in ['antibody_heavy', 'antibody_light']]

        total_antigen_residues = sum(
            chain_info[cid]['length'] for cid in antigen_chains if cid in chain_info
        )

        # Rule 1: No antigen chains
        if not antigen_chains or antigen_chain_str == 'NA' or not antigen_chain_str:
            return FilterResult(
                pdb_id=pdb_id,
                accepted=False,
                reason=FilterReason.NO_ANTIGEN,
                antigen_chains=antigen_chain_str,
                antigen_type=antigen_type_str,
                num_antigen_chains=0,
                total_antigen_residues=0,
                num_antibody_chains=len(antibody_chains),
                details="No antigen chains found or antigen_chain=NA"
            )

        # Rule 2: Check antigen type
        if antigen_type_str:
            antigen_types = [t.strip() for t in antigen_type_str.split('|')]

            # Reject hapten-only
            if all(t == 'hapten' for t in antigen_types if t):
                return FilterResult(
                    pdb_id=pdb_id,
                    accepted=False,
                    reason=FilterReason.HAPTEN_ONLY,
                    antigen_chains=antigen_chain_str,
                    antigen_type=antigen_type_str,
                    num_antigen_chains=len(antigen_chains),
                    total_antigen_residues=total_antigen_residues,
                    num_antibody_chains=len(antibody_chains),
                    details="Hapten-only antigen (not suitable for ESM-2)"
                )

            # Reject peptide-only (unless it's mixed with protein)
            has_protein = any(t == 'protein' for t in antigen_types)
            if not has_protein and any(t == 'peptide' for t in antigen_types):
                return FilterResult(
                    pdb_id=pdb_id,
                    accepted=False,
                    reason=FilterReason.PEPTIDE_ONLY,
                    antigen_chains=antigen_chain_str,
                    antigen_type=antigen_type_str,
                    num_antigen_chains=len(antigen_chains),
                    total_antigen_residues=total_antigen_residues,
                    num_antibody_chains=len(antibody_chains),
                    details="Peptide-only antigen (too short for reliable ESM-2 embedding)"
                )

            # Ensure at least one protein component
            if not has_protein:
                return FilterResult(
                    pdb_id=pdb_id,
                    accepted=False,
                    reason=FilterReason.NO_PROTEIN_ANTIGEN,
                    antigen_chains=antigen_chain_str,
                    antigen_type=antigen_type_str,
                    num_antigen_chains=len(antigen_chains),
                    total_antigen_residues=total_antigen_residues,
                    num_antibody_chains=len(antibody_chains),
                    details=f"No protein antigen component (type: {antigen_type_str})"
                )

        # Rule 3: Minimum residue count (at least 10 residues total)
        if total_antigen_residues < 10:
            return FilterResult(
                pdb_id=pdb_id,
                accepted=False,
                reason=FilterReason.INSUFFICIENT_ANTIGEN_RESIDUES,
                antigen_chains=antigen_chain_str,
                antigen_type=antigen_type_str,
                num_antigen_chains=len(antigen_chains),
                total_antigen_residues=total_antigen_residues,
                num_antibody_chains=len(antibody_chains),
                details=f"Only {total_antigen_residues} antigen residues (minimum: 10)"
            )

        # ACCEPT
        return FilterResult(
            pdb_id=pdb_id,
            accepted=True,
            reason=FilterReason.ACCEPTED,
            antigen_chains=antigen_chain_str,
            antigen_type=antigen_type_str,
            num_antigen_chains=len(antigen_chains),
            total_antigen_residues=total_antigen_residues,
            num_antibody_chains=len(antibody_chains),
            details=f"Protein antigen with {total_antigen_residues} residues"
        )

    def clean_structure(
        self, cif_path: Path, skip_filtering: bool = False
    ) -> Tuple[Optional[CleanedStructure], FilterResult]:
        """
        Clean a single CIF structure with filtering.

        Args:
            cif_path: Path to input mmCIF file
            skip_filtering: If True, skip antigen type filtering (for testing)

        Returns:
            Tuple of (CleanedStructure or None, FilterResult)
            - CleanedStructure is None if structure was filtered out
            - FilterResult contains filtering decision and details

        Raises:
            StructureCleaningError: If cleaning fails
        """
        if not cif_path.exists():
            raise StructureCleaningError(f"CIF file not found: {cif_path}")

        try:
            # Extract PDB ID from filename, removing common suffixes like '_cleaned'
            pdb_id = cif_path.stem.upper()
            for suffix in ['_CLEANED', '_CLEAN', '-CLEANED', '-CLEAN']:
                if pdb_id.endswith(suffix):
                    pdb_id = pdb_id[:-len(suffix)]
                    break

            # Step 1: Parse CIF with Biopython to get sequences
            chain_info = self._extract_chain_info(cif_path)

            # Step 2: Classify chains as antigen/antibody (using SAbDab metadata if available)
            chain_types = self._classify_chains(chain_info, pdb_id)

            # Step 3: Check if structure should be processed
            filter_result = self.should_process_structure(pdb_id, chain_info, chain_types)

            if not skip_filtering and not filter_result.accepted:
                # Structure filtered out - return None and filter result
                return None, filter_result

            # Step 4: Create standardized chain mappings
            chain_mappings = self._create_chain_mappings(
                pdb_id, chain_info, chain_types
            )

            # Step 5: Clean CIF with Gemmi and save
            output_path = self._clean_and_save_cif(
                cif_path, pdb_id, chain_mappings
            )

            cleaned = CleanedStructure(
                pdb_id=pdb_id,
                file_path=output_path,
                chain_mappings=chain_mappings,
                antigen_cluster_id=None  # Assigned later in pipeline
            )

            return cleaned, filter_result

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
                    'antigen_type': row.get('antigen_type', ''),
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


def save_filter_log(filter_results: List[FilterResult], output_path: Path):
    """
    Save filtering results to CSV file.

    Args:
        filter_results: List of FilterResult objects
        output_path: Path to output CSV file
    """
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)

        # Header
        writer.writerow([
            'pdb_id',
            'accepted',
            'filter_reason',
            'antigen_chains',
            'antigen_type',
            'num_antigen_chains',
            'total_antigen_residues',
            'num_antibody_chains',
            'details'
        ])

        # Data rows
        for result in filter_results:
            writer.writerow([
                result.pdb_id,
                'Yes' if result.accepted else 'No',
                result.reason.value,
                result.antigen_chains,
                result.antigen_type,
                result.num_antigen_chains,
                result.total_antigen_residues,
                result.num_antibody_chains,
                result.details
            ])


def generate_filter_summary(filter_results: List[FilterResult]) -> str:
    """
    Generate human-readable summary of filtering results.

    Args:
        filter_results: List of FilterResult objects

    Returns:
        Formatted summary string
    """
    total = len(filter_results)
    accepted = sum(1 for r in filter_results if r.accepted)
    rejected = total - accepted

    # Count by reason
    from collections import Counter
    reason_counts = Counter(r.reason for r in filter_results)

    summary = []
    summary.append("=" * 70)
    summary.append("Structure Filtering Summary")
    summary.append("=" * 70)
    summary.append(f"Total structures processed: {total}")
    summary.append(f"Accepted: {accepted} ({accepted/total*100:.1f}%)")
    summary.append(f"Rejected: {rejected} ({rejected/total*100:.1f}%)")
    summary.append("")
    summary.append("Rejection reasons:")

    for reason, count in reason_counts.most_common():
        if reason != FilterReason.ACCEPTED:
            percentage = count / total * 100
            summary.append(f"  - {reason.value}: {count} ({percentage:.1f}%)")

    summary.append("")
    summary.append(f"Accepted structures have protein antigens suitable for ESM-2 embedding")
    summary.append("=" * 70)

    return "\n".join(summary)
