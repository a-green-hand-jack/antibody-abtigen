"""
Epitope extraction module for the epitope-centric pipeline.

This module implements Phase 2: Epitope Definition.

Key tasks:
1. Extract epitope residues via distance-based contacts
2. Support multi-chain epitopes (e.g., trimers)
3. Merge epitope residues from multiple antigen chains
4. Calculate contact statistics

Reuses logic from:
- antibody_abtigen.epitope.get_epitope_residues()
- Adds multi-chain support and better metadata
"""

from pathlib import Path
from typing import Iterator, List, Tuple, Dict, Set

from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.Structure import Structure

from .core import (
    EpitopeExtractor,
    EpitopeResidues,
    CleanedStructure,
    ChainMapping,
    EpitopeExtractionError,
)


class GeometricEpitopeExtractor(EpitopeExtractor):
    """
    Geometric epitope extractor using distance-based contacts.

    An epitope is defined as the set of antigen residues that have
    at least one atom within a distance threshold of any antibody atom.

    Supports:
    - Multi-chain epitopes (antibody binding at chain interfaces)
    - Customizable distance thresholds
    - Contact statistics

    Usage:
        >>> extractor = GeometricEpitopeExtractor(distance_threshold=5.0)
        >>> epitope = extractor.extract_epitope(cleaned_structure)
        >>> print(epitope.total_residue_count())
        23
    """

    def __init__(self, distance_threshold: float = 5.0):
        """
        Initialize extractor.

        Args:
            distance_threshold: Contact distance cutoff in Angstroms (default: 5.0)
        """
        self.distance_threshold = distance_threshold

    def extract_epitope(
        self, cleaned_structure: CleanedStructure
    ) -> EpitopeResidues:
        """
        Extract epitope from a cleaned structure.

        Args:
            cleaned_structure: Cleaned structure with standardized chains

        Returns:
            EpitopeResidues with multi-chain support

        Raises:
            EpitopeExtractionError: If extraction fails
        """
        try:
            # Parse structure
            structure = self._parse_structure(cleaned_structure.file_path)

            # Get antigen and antibody chains
            antigen_chains = [
                m.standardized_chain_id
                for m in cleaned_structure.get_antigen_chains()
            ]
            antibody_chains = [
                m.standardized_chain_id
                for m in cleaned_structure.get_antibody_chains()
            ]

            if not antigen_chains:
                raise EpitopeExtractionError(
                    f"No antigen chains in {cleaned_structure.pdb_id}"
                )

            if not antibody_chains:
                raise EpitopeExtractionError(
                    f"No antibody chains in {cleaned_structure.pdb_id}"
                )

            # Extract epitope residues
            epitope_residues_by_chain = self._get_epitope_residues(
                structure, antigen_chains, antibody_chains
            )

            # Count total contacts
            total_contacts = sum(
                len(residues) for residues in epitope_residues_by_chain.values()
            )

            # Create epitope ID
            epitope_id = f"{cleaned_structure.pdb_id.lower()}_epi"

            return EpitopeResidues(
                epitope_id=epitope_id,
                pdb_id=cleaned_structure.pdb_id,
                antigen_chains=epitope_residues_by_chain,
                antibody_chains=antibody_chains,
                distance_threshold=self.distance_threshold,
                num_contacts=total_contacts
            )

        except Exception as e:
            raise EpitopeExtractionError(
                f"Failed to extract epitope from {cleaned_structure.pdb_id}: {str(e)}"
            ) from e

    def batch_extract(
        self, cleaned_structures: List[CleanedStructure]
    ) -> Iterator[Tuple[CleanedStructure, EpitopeResidues]]:
        """
        Extract epitopes from multiple structures.

        Args:
            cleaned_structures: List of cleaned structures

        Yields:
            Tuples of (CleanedStructure, EpitopeResidues)
        """
        for structure in cleaned_structures:
            try:
                epitope = self.extract_epitope(structure)
                yield structure, epitope
            except EpitopeExtractionError as e:
                # Log error but continue
                print(f"Warning: {e}")
                continue

    def _parse_structure(self, file_path: Path) -> Structure:
        """
        Parse a structure file (CIF or PDB).

        Args:
            file_path: Path to structure file

        Returns:
            Biopython Structure object

        Raises:
            EpitopeExtractionError: If parsing fails
        """
        if not file_path.exists():
            raise EpitopeExtractionError(f"File not found: {file_path}")

        try:
            if file_path.suffix == '.cif':
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)

            structure = parser.get_structure("structure", str(file_path))
            return structure

        except Exception as e:
            raise EpitopeExtractionError(
                f"Failed to parse {file_path}: {str(e)}"
            ) from e

    def _get_epitope_residues(
        self,
        structure: Structure,
        antigen_chain_ids: List[str],
        antibody_chain_ids: List[str]
    ) -> Dict[str, List[int]]:
        """
        Extract epitope residues from structure.

        Uses NeighborSearch for efficient spatial queries.

        Args:
            structure: Biopython Structure
            antigen_chain_ids: Standardized antigen chain IDs
            antibody_chain_ids: Standardized antibody chain IDs

        Returns:
            Dict mapping chain_id → [auth_seq_id, ...]
        """
        model = structure[0]

        # Collect all antibody atoms
        antibody_atoms = []
        for chain_id in antibody_chain_ids:
            if chain_id not in model:
                continue
            chain = model[chain_id]
            for residue in chain:
                # Only standard residues
                if residue.get_id()[0] == ' ':
                    for atom in residue:
                        antibody_atoms.append(atom)

        if not antibody_atoms:
            raise EpitopeExtractionError("No antibody atoms found")

        # Build neighbor search
        ns = NeighborSearch(antibody_atoms)

        # Find antigen residues in contact
        epitope_by_chain: Dict[str, Set[int]] = {}

        for chain_id in antigen_chain_ids:
            if chain_id not in model:
                continue

            chain = model[chain_id]
            epitope_residues = set()

            for residue in chain:
                # Only standard residues
                if residue.get_id()[0] != ' ':
                    continue

                # Check if any atom is near antibody
                for atom in residue:
                    nearby = ns.search(
                        atom.get_coord(),
                        self.distance_threshold
                    )
                    if nearby:
                        # Use auth_seq_id (PDB numbering)
                        epitope_residues.add(residue.get_id()[1])
                        break  # Found contact, move to next residue

            if epitope_residues:
                epitope_by_chain[chain_id] = epitope_residues

        # Convert sets to sorted lists
        return {
            chain_id: sorted(list(residues))
            for chain_id, residues in epitope_by_chain.items()
        }
