"""
Structure-based epitope validation and pocket expansion module.

This module implements:
1. Pocket Cropper: Expand epitope residues to contiguous pocket regions
2. Structure Validator: Validate groups by structural RMSD of pocket regions

Inspired by Boltz-2 Algorithm 3: Affinity Cropper

Key concepts:
- Epitope: Residues in direct contact with antibody (< 5Å)
- Pocket: Expanded region including epitope + neighboring residues
- Structure similarity: Based on pocket Cα RMSD after optimal superposition
"""

import json
import logging
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set, Any
import csv

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Structure import Structure

from .core import (
    CleanedStructure,
    EpitopeResidues,
)
from .grouper import GroupingOutput, GroupResult, GroupMember

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class PocketResidues:
    """
    Expanded pocket region around epitope.

    Unlike EpitopeResidues which only contains contact residues,
    PocketResidues includes neighboring residues for more complete
    structural comparison.
    """
    epitope_id: str
    pdb_id: str
    # Pocket residues by chain: chain_id -> sorted list of auth_seq_ids
    pocket_chains: Dict[str, List[int]]
    # Original epitope residues (subset of pocket)
    epitope_chains: Dict[str, List[int]]
    # Parameters used
    neighborhood_size: int
    max_pocket_residues: int

    def total_pocket_residues(self) -> int:
        """Total number of residues in pocket."""
        return sum(len(res) for res in self.pocket_chains.values())

    def total_epitope_residues(self) -> int:
        """Total number of original epitope residues."""
        return sum(len(res) for res in self.epitope_chains.values())

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PocketCoordinates:
    """3D coordinates of pocket Cα atoms."""
    epitope_id: str
    pdb_id: str
    # Ordered list of (chain_id, res_id, coord) for each Cα
    coords: List[Tuple[str, int, np.ndarray]]

    def get_coords_array(self) -> np.ndarray:
        """Get Cα coordinates as (N, 3) array."""
        if not self.coords:
            return np.empty((0, 3))
        return np.array([c[2] for c in self.coords])

    def get_residue_keys(self) -> List[Tuple[str, int]]:
        """Get (chain_id, res_id) keys for matching."""
        return [(c[0], c[1]) for c in self.coords]


@dataclass
class StructuralValidationResult:
    """Result of validating a single group member."""
    epitope_id: str
    pdb_id: str
    rmsd_to_reference: Optional[float]  # None for reference
    atoms_aligned: int
    atoms_in_pocket: int
    alignment_coverage: float  # atoms_aligned / atoms_in_pocket
    passes_threshold: bool
    is_reference: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ValidatedGroupResult:
    """A validated epitope group with structural metrics."""
    group_id: str
    original_member_count: int
    validated_member_count: int
    reference_epitope_id: str
    reference_pdb_id: str
    avg_rmsd: float
    max_rmsd: float
    min_coverage: float
    avg_coverage: float
    validated_members: List[StructuralValidationResult]
    rejected_members: List[StructuralValidationResult]
    # Original embedding similarity (for comparison)
    embedding_similarity: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        result = {
            'group_id': self.group_id,
            'original_member_count': self.original_member_count,
            'validated_member_count': self.validated_member_count,
            'reference_epitope_id': self.reference_epitope_id,
            'reference_pdb_id': self.reference_pdb_id,
            'avg_rmsd': self.avg_rmsd,
            'max_rmsd': self.max_rmsd,
            'min_coverage': self.min_coverage,
            'avg_coverage': self.avg_coverage,
            'embedding_similarity': self.embedding_similarity,
            'validated_members': [m.to_dict() for m in self.validated_members],
            'rejected_members': [m.to_dict() for m in self.rejected_members],
        }
        return result


@dataclass
class ValidationOutput:
    """Complete output from the structure validator."""
    validated_groups: List[ValidatedGroupResult]
    total_groups_input: int
    total_groups_output: int
    groups_unchanged: int
    groups_reduced: int
    groups_eliminated: int
    rmsd_threshold: float
    min_coverage: float
    min_group_size: int
    neighborhood_size: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            'total_groups_input': self.total_groups_input,
            'total_groups_output': self.total_groups_output,
            'groups_unchanged': self.groups_unchanged,
            'groups_reduced': self.groups_reduced,
            'groups_eliminated': self.groups_eliminated,
            'rmsd_threshold': self.rmsd_threshold,
            'min_coverage': self.min_coverage,
            'min_group_size': self.min_group_size,
            'neighborhood_size': self.neighborhood_size,
            'validated_groups': [g.to_dict() for g in self.validated_groups],
        }


# =============================================================================
# Pocket Cropper (inspired by Boltz-2 Algorithm 3)
# =============================================================================

class PocketCropper:
    """
    Expand epitope residues to contiguous pocket regions.

    Algorithm (adapted from Boltz-2 Affinity Cropper):
    1. Start with epitope contact residues
    2. For each contact residue, expand to include neighbors
    3. Ensure contiguous regions of at least neighborhood_size residues
    4. Cap at max_pocket_residues to prevent huge pockets

    This creates a more complete binding site representation for
    structure-based comparison.
    """

    def __init__(
        self,
        neighborhood_size: int = 10,
        max_pocket_residues: int = 50,
    ):
        """
        Initialize pocket cropper.

        Args:
            neighborhood_size: Minimum window size around each contact residue
            max_pocket_residues: Maximum total pocket residues per chain
        """
        self.neighborhood_size = neighborhood_size
        self.max_pocket_residues = max_pocket_residues

    def expand_epitope_to_pocket(
        self,
        epitope: EpitopeResidues,
        structure: Structure,
    ) -> PocketResidues:
        """
        Expand epitope residues to pocket regions.

        Args:
            epitope: Original epitope with contact residues
            structure: Biopython Structure object

        Returns:
            PocketResidues with expanded pocket regions
        """
        model = structure[0]
        pocket_chains: Dict[str, List[int]] = {}

        for chain_id, contact_residues in epitope.antigen_chains.items():
            if chain_id not in model:
                continue

            chain = model[chain_id]

            # Get all residue IDs in this chain (sorted)
            all_res_ids = sorted([
                res.get_id()[1]
                for res in chain
                if res.get_id()[0] == ' '  # Standard residues only
            ])

            if not all_res_ids:
                continue

            # Expand each contact residue to a window
            pocket_residues: Set[int] = set()

            # Sort contact residues by position
            sorted_contacts = sorted(contact_residues)

            for contact_res in sorted_contacts:
                if contact_res not in all_res_ids:
                    continue

                # Find index of contact residue
                try:
                    center_idx = all_res_ids.index(contact_res)
                except ValueError:
                    continue

                # Expand window until we have neighborhood_size residues
                half_window = self.neighborhood_size // 2
                min_idx = max(0, center_idx - half_window)
                max_idx = min(len(all_res_ids), center_idx + half_window + 1)

                # Add residues in window
                window_residues = all_res_ids[min_idx:max_idx]
                pocket_residues.update(window_residues)

                # Check if we've hit the limit
                if len(pocket_residues) >= self.max_pocket_residues:
                    break

            # Sort and limit
            pocket_list = sorted(pocket_residues)[:self.max_pocket_residues]

            if pocket_list:
                pocket_chains[chain_id] = pocket_list

        return PocketResidues(
            epitope_id=epitope.epitope_id,
            pdb_id=epitope.pdb_id,
            pocket_chains=pocket_chains,
            epitope_chains=epitope.antigen_chains,
            neighborhood_size=self.neighborhood_size,
            max_pocket_residues=self.max_pocket_residues,
        )

    def extract_pocket_coordinates(
        self,
        pocket: PocketResidues,
        structure: Structure,
    ) -> PocketCoordinates:
        """
        Extract Cα coordinates for pocket residues.

        Args:
            pocket: Pocket residues definition
            structure: Biopython Structure object

        Returns:
            PocketCoordinates with Cα positions
        """
        model = structure[0]
        coords: List[Tuple[str, int, np.ndarray]] = []

        for chain_id in sorted(pocket.pocket_chains.keys()):
            res_ids = pocket.pocket_chains[chain_id]

            if chain_id not in model:
                continue

            chain = model[chain_id]

            for res_id in res_ids:
                # Find residue
                for residue in chain:
                    if residue.get_id()[1] == res_id and residue.get_id()[0] == ' ':
                        if 'CA' in residue:
                            ca_coord = residue['CA'].get_coord()
                            coords.append((chain_id, res_id, ca_coord))
                        break

        return PocketCoordinates(
            epitope_id=pocket.epitope_id,
            pdb_id=pocket.pdb_id,
            coords=coords,
        )


# =============================================================================
# Structure Validator
# =============================================================================

class StructureValidator:
    """
    Validate epitope groups based on structural alignment of pocket regions.

    Strategy:
    1. Expand epitopes to pocket regions (using PocketCropper)
    2. For each group, compute pairwise pocket Cα RMSD
    3. Filter members exceeding RMSD threshold
    4. Require minimum alignment coverage
    5. Optionally re-select reference as most structurally central member
    """

    def __init__(
        self,
        rmsd_threshold: float = 3.0,
        min_coverage: float = 0.5,
        min_group_size: int = 2,
        neighborhood_size: int = 10,
        max_pocket_residues: int = 50,
        reselect_reference: bool = True,
    ):
        """
        Initialize structure validator.

        Args:
            rmsd_threshold: Maximum RMSD in Angstroms for structural similarity
            min_coverage: Minimum fraction of pocket atoms that must align
            min_group_size: Minimum members for a valid group after filtering
            neighborhood_size: Window size for pocket expansion
            max_pocket_residues: Maximum pocket size per chain
            reselect_reference: Whether to re-select reference based on structure
        """
        self.rmsd_threshold = rmsd_threshold
        self.min_coverage = min_coverage
        self.min_group_size = min_group_size
        self.neighborhood_size = neighborhood_size
        self.max_pocket_residues = max_pocket_residues
        self.reselect_reference = reselect_reference

        self.cropper = PocketCropper(
            neighborhood_size=neighborhood_size,
            max_pocket_residues=max_pocket_residues,
        )

    def validate_groups(
        self,
        grouping_output: GroupingOutput,
        epitopes: Dict[str, EpitopeResidues],
        structures: Dict[str, CleanedStructure],
    ) -> ValidationOutput:
        """
        Validate all groups using structural alignment.

        Args:
            grouping_output: Output from embedding-based grouping
            epitopes: Dict of epitope_id -> EpitopeResidues
            structures: Dict of pdb_id -> CleanedStructure

        Returns:
            ValidationOutput with validated groups
        """
        validated_groups: List[ValidatedGroupResult] = []
        groups_unchanged = 0
        groups_reduced = 0
        groups_eliminated = 0

        # Parse structures and compute pocket coordinates
        logger.info("Expanding epitopes to pocket regions...")
        pocket_coords = self._compute_all_pocket_coords(epitopes, structures)

        logger.info(f"Validating {len(grouping_output.groups)} groups...")

        for group in grouping_output.groups:
            validated = self._validate_single_group(group, pocket_coords)

            if validated is None:
                groups_eliminated += 1
            elif validated.validated_member_count == group.member_count:
                groups_unchanged += 1
                validated_groups.append(validated)
            else:
                groups_reduced += 1
                validated_groups.append(validated)

        return ValidationOutput(
            validated_groups=validated_groups,
            total_groups_input=len(grouping_output.groups),
            total_groups_output=len(validated_groups),
            groups_unchanged=groups_unchanged,
            groups_reduced=groups_reduced,
            groups_eliminated=groups_eliminated,
            rmsd_threshold=self.rmsd_threshold,
            min_coverage=self.min_coverage,
            min_group_size=self.min_group_size,
            neighborhood_size=self.neighborhood_size,
        )

    def _compute_all_pocket_coords(
        self,
        epitopes: Dict[str, EpitopeResidues],
        structures: Dict[str, CleanedStructure],
    ) -> Dict[str, PocketCoordinates]:
        """Compute pocket coordinates for all epitopes."""
        pocket_coords: Dict[str, PocketCoordinates] = {}

        # Cache parsed structures
        parsed_structures: Dict[str, Structure] = {}

        for epitope_id, epitope in epitopes.items():
            pdb_id = epitope.pdb_id

            # Get or parse structure
            if pdb_id not in parsed_structures:
                if pdb_id not in structures:
                    logger.warning(f"Structure not found for {pdb_id}")
                    continue

                cleaned = structures[pdb_id]
                parsed_structures[pdb_id] = self._parse_structure(cleaned.file_path)

            structure = parsed_structures[pdb_id]

            # Expand to pocket
            pocket = self.cropper.expand_epitope_to_pocket(epitope, structure)

            # Extract coordinates
            coords = self.cropper.extract_pocket_coordinates(pocket, structure)

            if coords.coords:
                pocket_coords[epitope_id] = coords

        logger.info(f"Computed pocket coordinates for {len(pocket_coords)} epitopes")
        return pocket_coords

    def _validate_single_group(
        self,
        group: GroupResult,
        pocket_coords: Dict[str, PocketCoordinates],
    ) -> Optional[ValidatedGroupResult]:
        """
        Validate a single group using structural alignment.

        Returns None if group is eliminated (< min_group_size valid members).
        """
        # Get reference coordinates
        ref_epitope_id = group.reference_epitope_id
        if ref_epitope_id not in pocket_coords:
            logger.warning(f"Reference {ref_epitope_id} not found in pocket coords")
            return None

        ref_coords = pocket_coords[ref_epitope_id]
        ref_array = ref_coords.get_coords_array()

        if len(ref_array) < 3:
            logger.warning(f"Reference {ref_epitope_id} has too few atoms ({len(ref_array)})")
            return None

        validated_members: List[StructuralValidationResult] = []
        rejected_members: List[StructuralValidationResult] = []

        # Add reference as validated
        ref_result = StructuralValidationResult(
            epitope_id=ref_epitope_id,
            pdb_id=ref_coords.pdb_id,
            rmsd_to_reference=None,
            atoms_aligned=len(ref_array),
            atoms_in_pocket=len(ref_array),
            alignment_coverage=1.0,
            passes_threshold=True,
            is_reference=True,
        )
        validated_members.append(ref_result)

        # Validate each member
        for member in group.members:
            if member.is_reference:
                continue

            epitope_id = member.epitope_id

            if epitope_id not in pocket_coords:
                # No coordinates - reject
                rejected_members.append(StructuralValidationResult(
                    epitope_id=epitope_id,
                    pdb_id=member.pdb_id,
                    rmsd_to_reference=None,
                    atoms_aligned=0,
                    atoms_in_pocket=0,
                    alignment_coverage=0.0,
                    passes_threshold=False,
                    is_reference=False,
                ))
                continue

            mob_coords = pocket_coords[epitope_id]
            mob_array = mob_coords.get_coords_array()

            if len(mob_array) < 3:
                rejected_members.append(StructuralValidationResult(
                    epitope_id=epitope_id,
                    pdb_id=member.pdb_id,
                    rmsd_to_reference=None,
                    atoms_aligned=len(mob_array),
                    atoms_in_pocket=len(mob_array),
                    alignment_coverage=0.0,
                    passes_threshold=False,
                    is_reference=False,
                ))
                continue

            # Compute RMSD
            rmsd, n_aligned = self._compute_rmsd(ref_array, mob_array)

            # Calculate coverage (fraction of smaller pocket that aligned)
            min_pocket_size = min(len(ref_array), len(mob_array))
            coverage = n_aligned / min_pocket_size if min_pocket_size > 0 else 0.0

            # Check thresholds
            passes = (
                rmsd <= self.rmsd_threshold and
                coverage >= self.min_coverage
            )

            result = StructuralValidationResult(
                epitope_id=epitope_id,
                pdb_id=member.pdb_id,
                rmsd_to_reference=rmsd,
                atoms_aligned=n_aligned,
                atoms_in_pocket=len(mob_array),
                alignment_coverage=coverage,
                passes_threshold=passes,
                is_reference=False,
            )

            if passes:
                validated_members.append(result)
            else:
                rejected_members.append(result)

        # Check minimum group size
        if len(validated_members) < self.min_group_size:
            return None

        # Calculate statistics
        valid_rmsds = [
            m.rmsd_to_reference
            for m in validated_members
            if m.rmsd_to_reference is not None
        ]
        valid_coverages = [m.alignment_coverage for m in validated_members]

        avg_rmsd = np.mean(valid_rmsds) if valid_rmsds else 0.0
        max_rmsd = max(valid_rmsds) if valid_rmsds else 0.0
        min_coverage = min(valid_coverages) if valid_coverages else 0.0
        avg_coverage = np.mean(valid_coverages) if valid_coverages else 0.0

        return ValidatedGroupResult(
            group_id=group.group_id,
            original_member_count=group.member_count,
            validated_member_count=len(validated_members),
            reference_epitope_id=ref_epitope_id,
            reference_pdb_id=ref_coords.pdb_id,
            avg_rmsd=float(avg_rmsd),
            max_rmsd=float(max_rmsd),
            min_coverage=float(min_coverage),
            avg_coverage=float(avg_coverage),
            validated_members=validated_members,
            rejected_members=rejected_members,
            embedding_similarity=group.avg_similarity,
        )

    def _compute_rmsd(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray,
    ) -> Tuple[float, int]:
        """
        Compute RMSD between two coordinate sets using SVD superposition.

        Uses Kabsch algorithm for optimal rotation.
        Handles different sizes by using iterative closest point matching.

        Args:
            coords1: Reference coordinates (N, 3)
            coords2: Mobile coordinates (M, 3)

        Returns:
            Tuple of (RMSD, number of atoms aligned)
        """
        n1 = len(coords1)
        n2 = len(coords2)

        if n1 == 0 or n2 == 0:
            return float('inf'), 0

        # Use the smaller set size for alignment
        n_align = min(n1, n2)

        if n_align < 3:
            return float('inf'), n_align

        # If sizes differ, match closest points
        if n1 != n2:
            # Use subset of larger to match smaller
            if n1 > n2:
                # Find closest points in coords1 to each point in coords2
                matched_indices = self._match_closest_points(coords1, coords2)
                fixed = coords1[matched_indices]
                moving = coords2
            else:
                matched_indices = self._match_closest_points(coords2, coords1)
                fixed = coords1
                moving = coords2[matched_indices]
        else:
            fixed = coords1
            moving = coords2

        # Center coordinates
        fixed_center = fixed.mean(axis=0)
        moving_center = moving.mean(axis=0)
        fixed_centered = fixed - fixed_center
        moving_centered = moving - moving_center

        # SVD for optimal rotation (Kabsch algorithm)
        H = moving_centered.T @ fixed_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Handle reflection
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # Apply rotation and compute RMSD
        moving_aligned = moving_centered @ R
        diff = fixed_centered - moving_aligned
        rmsd = np.sqrt((diff ** 2).sum() / len(diff))

        return float(rmsd), len(fixed)

    def _match_closest_points(
        self,
        large_coords: np.ndarray,
        small_coords: np.ndarray,
    ) -> np.ndarray:
        """
        Match points in larger set to closest points in smaller set.

        Returns indices into large_coords that best match small_coords.
        """
        from scipy.spatial import cKDTree

        tree = cKDTree(large_coords)
        distances, indices = tree.query(small_coords, k=1)

        # Remove duplicates (each point in large can only be matched once)
        used = set()
        unique_indices = []
        for idx in indices:
            if idx not in used:
                unique_indices.append(idx)
                used.add(idx)
            else:
                # Find next closest unused point
                for i in range(len(large_coords)):
                    if i not in used:
                        unique_indices.append(i)
                        used.add(i)
                        break

        return np.array(unique_indices[:len(small_coords)])

    def _parse_structure(self, file_path: Path) -> Structure:
        """Parse a structure file."""
        if file_path.suffix == '.cif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        return parser.get_structure("structure", str(file_path))


# =============================================================================
# I/O Functions
# =============================================================================

def save_validated_groups_json(
    output: ValidationOutput,
    output_path: Path,
) -> None:
    """Save validated groups to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output.to_dict(), f, indent=2, default=str)
    logger.info(f"Saved validated groups to {output_path}")


def save_validation_report_csv(
    output: ValidationOutput,
    output_path: Path,
) -> None:
    """Save validation report as CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for group in output.validated_groups:
        for member in group.validated_members + group.rejected_members:
            rows.append({
                'group_id': group.group_id,
                'epitope_id': member.epitope_id,
                'pdb_id': member.pdb_id,
                'is_reference': member.is_reference,
                'rmsd_to_reference': member.rmsd_to_reference,
                'atoms_aligned': member.atoms_aligned,
                'atoms_in_pocket': member.atoms_in_pocket,
                'alignment_coverage': member.alignment_coverage,
                'passes_threshold': member.passes_threshold,
                'rmsd_threshold': output.rmsd_threshold,
                'min_coverage_threshold': output.min_coverage,
            })

    if rows:
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        logger.info(f"Saved validation report to {output_path}")


def generate_validation_report(output: ValidationOutput) -> str:
    """Generate human-readable validation report."""
    lines = [
        "=" * 70,
        "Structure Validation Report",
        "=" * 70,
        "",
        f"Parameters:",
        f"  RMSD threshold:      {output.rmsd_threshold:.1f} Å",
        f"  Min coverage:        {output.min_coverage:.0%}",
        f"  Min group size:      {output.min_group_size}",
        f"  Neighborhood size:   {output.neighborhood_size} residues",
        "",
        f"Results:",
        f"  Groups input:        {output.total_groups_input}",
        f"  Groups output:       {output.total_groups_output}",
        f"  Groups unchanged:    {output.groups_unchanged}",
        f"  Groups reduced:      {output.groups_reduced}",
        f"  Groups eliminated:   {output.groups_eliminated}",
        "",
    ]

    if output.validated_groups:
        lines.append("Validated Groups:")
        lines.append("-" * 70)

        for group in output.validated_groups:
            lines.append(
                f"  {group.group_id}: {group.validated_member_count}/{group.original_member_count} members, "
                f"RMSD={group.avg_rmsd:.2f}Å (max={group.max_rmsd:.2f}Å), "
                f"coverage={group.avg_coverage:.0%}"
            )

    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)


def load_validated_groups_json(input_path: Path) -> ValidationOutput:
    """Load validated groups from JSON file."""
    with open(input_path, 'r') as f:
        data = json.load(f)

    validated_groups = []
    for g in data.get('validated_groups', []):
        validated_members = [
            StructuralValidationResult(**m)
            for m in g.get('validated_members', [])
        ]
        rejected_members = [
            StructuralValidationResult(**m)
            for m in g.get('rejected_members', [])
        ]

        validated_groups.append(ValidatedGroupResult(
            group_id=g['group_id'],
            original_member_count=g['original_member_count'],
            validated_member_count=g['validated_member_count'],
            reference_epitope_id=g['reference_epitope_id'],
            reference_pdb_id=g['reference_pdb_id'],
            avg_rmsd=g['avg_rmsd'],
            max_rmsd=g['max_rmsd'],
            min_coverage=g['min_coverage'],
            avg_coverage=g['avg_coverage'],
            validated_members=validated_members,
            rejected_members=rejected_members,
            embedding_similarity=g.get('embedding_similarity', 0.0),
        ))

    return ValidationOutput(
        validated_groups=validated_groups,
        total_groups_input=data.get('total_groups_input', 0),
        total_groups_output=data.get('total_groups_output', 0),
        groups_unchanged=data.get('groups_unchanged', 0),
        groups_reduced=data.get('groups_reduced', 0),
        groups_eliminated=data.get('groups_eliminated', 0),
        rmsd_threshold=data.get('rmsd_threshold', 3.0),
        min_coverage=data.get('min_coverage', 0.5),
        min_group_size=data.get('min_group_size', 2),
        neighborhood_size=data.get('neighborhood_size', 10),
    )
