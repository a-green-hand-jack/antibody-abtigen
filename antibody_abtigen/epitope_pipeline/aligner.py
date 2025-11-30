"""
Structure aligner module for the epitope-centric pipeline.

This module implements Phase 5: Epitope-based Structure Alignment.

Key features:
1. Align structures based on epitope Cα atoms only (not full chain)
2. Apply transformation to entire complex (antigen + antibody)
3. Save aligned structures as separate CIF files
4. Generate alignment metadata (RMSD, transformation matrices)

The alignment strategy:
- Reference structure stays fixed
- Mobile structures are aligned to reference based on epitope residues
- Antibody chains move with their antigen (preserving relative position)
"""

import json
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field, asdict

import numpy as np
import gemmi

from .core import (
    StructureAligner,
    EpitopeGroup,
    EpitopeResidues,
    CleanedStructure,
    AlignedComplex,
    AlignmentError,
)
from .grouper import GroupResult, GroupMember

logger = logging.getLogger(__name__)


# PyMOL detection (reuse logic from structure.py)
PYMOL_PYTHON_PATH: Optional[str] = None
PYMOL_AVAILABLE: bool = False


def _detect_pymol() -> Tuple[bool, Optional[str]]:
    """
    Detect PyMOL installation.

    Returns:
        (available, python_path) tuple
    """
    import shutil

    # Check common conda environments
    conda_envs = [
        "/ibex/user/wuj0c/env/mamba/pymol-env/bin/python",
        "pymol-env",
        "pymol",
        "pymol-open-source",
    ]

    # Try conda envs
    for env_name in conda_envs:
        if env_name.startswith("/") and Path(env_name).exists():
            return True, env_name

    # Try mamba/conda run
    for env in ["pymol-env", "pymol"]:
        try:
            result = subprocess.run(
                ["mamba", "run", "-n", env, "python", "-c", "import pymol"],
                capture_output=True,
                timeout=30
            )
            if result.returncode == 0:
                # Get the python path
                result2 = subprocess.run(
                    ["mamba", "run", "-n", env, "which", "python"],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result2.returncode == 0:
                    return True, result2.stdout.strip()
        except (subprocess.TimeoutExpired, FileNotFoundError):
            continue

    # Try direct import
    try:
        import pymol
        return True, None  # Use direct import
    except ImportError:
        pass

    return False, None


def _init_pymol():
    """Initialize PyMOL detection."""
    global PYMOL_AVAILABLE, PYMOL_PYTHON_PATH
    if PYMOL_AVAILABLE:
        return
    PYMOL_AVAILABLE, PYMOL_PYTHON_PATH = _detect_pymol()
    if PYMOL_AVAILABLE:
        logger.info(f"PyMOL available: python_path={PYMOL_PYTHON_PATH or 'direct import'}")
    else:
        logger.warning("PyMOL not available, will use Biopython fallback")


@dataclass
class AlignmentResult:
    """Result of aligning a single structure to reference."""
    epitope_id: str
    pdb_id: str
    rmsd: float
    num_atoms_aligned: int
    rotation_matrix: List[List[float]]  # 3x3
    translation_vector: List[float]  # 3
    is_reference: bool = False


@dataclass
class GroupAlignmentOutput:
    """Complete output from aligning a group."""
    group_id: str
    reference_pdb: str
    reference_epitope_id: str
    members: List[AlignmentResult]
    avg_rmsd: float
    output_dir: Path
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'group_id': self.group_id,
            'reference_pdb': self.reference_pdb,
            'reference_epitope_id': self.reference_epitope_id,
            'avg_rmsd': self.avg_rmsd,
            'members': [
                {
                    'pdb_id': m.pdb_id,
                    'epitope_id': m.epitope_id,
                    'is_reference': m.is_reference,
                    'rmsd': m.rmsd if not m.is_reference else None,
                    'num_atoms_aligned': m.num_atoms_aligned,
                }
                for m in self.members
            ],
            'output_dir': str(self.output_dir),
        }


class PyMOLStructureAligner(StructureAligner):
    """
    PyMOL-based structure aligner for epitope-centric alignment.

    Aligns structures based on epitope Cα atoms using PyMOL's super command.
    Falls back to Biopython Superimposer if PyMOL is unavailable.

    Usage:
        >>> aligner = PyMOLStructureAligner()
        >>> results = aligner.align_group(group, epitopes, structures, output_dir)
    """

    def __init__(
        self,
        use_super: bool = True,
        fallback_to_biopython: bool = True
    ):
        """
        Initialize aligner.

        Args:
            use_super: Use PyMOL super (True) or align (False)
            fallback_to_biopython: Fall back to Biopython if PyMOL unavailable
        """
        self.use_super = use_super
        self.fallback_to_biopython = fallback_to_biopython
        _init_pymol()

    def align_group(
        self,
        group: EpitopeGroup,
        epitopes: Dict[str, EpitopeResidues],
        structures: Dict[str, CleanedStructure],
        output_dir: Path
    ) -> List[AlignedComplex]:
        """
        Align all structures in a group to the reference.

        This is the interface method required by StructureAligner.

        Args:
            group: EpitopeGroup with reference and members
            epitopes: Mapping epitope_id -> EpitopeResidues
            structures: Mapping pdb_id -> CleanedStructure
            output_dir: Where to save aligned CIF files

        Returns:
            List of AlignedComplex objects
        """
        # Convert to GroupResult format for internal processing
        group_result = GroupResult(
            group_id=group.group_id,
            reference_epitope_id=group.reference_epitope_id,
            member_count=len(group.member_epitope_ids),
            avg_similarity=group.avg_similarity,
            min_similarity=0.0,
            max_similarity=1.0,
            members=[
                GroupMember(
                    epitope_id=eid,
                    pdb_id=eid.split('_')[0].upper(),  # Extract PDB ID from epitope_id
                    is_reference=(eid == group.reference_epitope_id),
                    antigen_chains=[],
                    epitope_residues={},
                    total_residues=0
                )
                for eid in group.member_epitope_ids
            ]
        )

        output = self.align_group_detailed(
            group_result, epitopes, structures, output_dir
        )

        # Convert to AlignedComplex format
        aligned_complexes = []
        for member in output.members:
            if member.is_reference:
                # Reference files
                antigen_path = output.output_dir / "reference" / f"{member.pdb_id}_antigen.cif"
                antibody_path = output.output_dir / "reference" / f"{member.pdb_id}_antibody.cif"
            else:
                # Aligned files
                antigen_path = output.output_dir / "aligned" / f"{member.pdb_id}_antigen.cif"
                antibody_path = output.output_dir / "aligned" / f"{member.pdb_id}_antibody.cif"

            # Create AlignedComplex for antigen
            aligned_complexes.append(AlignedComplex(
                epitope_id=member.epitope_id,
                pdb_id=member.pdb_id,
                group_id=output.group_id,
                aligned_file_path=antigen_path,
                rmsd_to_reference=member.rmsd if not member.is_reference else 0.0,
                alignment_method='pymol_super' if PYMOL_AVAILABLE else 'biopython',
                is_reference=member.is_reference,
                metadata={'type': 'antigen'}
            ))

        return aligned_complexes

    def align_group_detailed(
        self,
        group: GroupResult,
        epitopes: Dict[str, EpitopeResidues],
        structures: Dict[str, CleanedStructure],
        output_dir: Path
    ) -> GroupAlignmentOutput:
        """
        Align all structures in a group with detailed output.

        Args:
            group: GroupResult from grouper
            epitopes: Mapping epitope_id -> EpitopeResidues
            structures: Mapping pdb_id -> CleanedStructure
            output_dir: Where to save aligned CIF files

        Returns:
            GroupAlignmentOutput with all alignment results
        """
        group_dir = output_dir / group.group_id
        ref_dir = group_dir / "reference"
        aligned_dir = group_dir / "aligned"

        ref_dir.mkdir(parents=True, exist_ok=True)
        aligned_dir.mkdir(parents=True, exist_ok=True)

        # Find reference
        reference_member = None
        for member in group.members:
            if member.is_reference:
                reference_member = member
                break

        if reference_member is None:
            raise AlignmentError(f"No reference found in group {group.group_id}")

        ref_epitope_id = reference_member.epitope_id
        ref_pdb_id = reference_member.pdb_id

        if ref_pdb_id not in structures:
            raise AlignmentError(f"Reference structure {ref_pdb_id} not found")
        if ref_epitope_id not in epitopes:
            raise AlignmentError(f"Reference epitope {ref_epitope_id} not found")

        ref_structure = structures[ref_pdb_id]
        ref_epitope = epitopes[ref_epitope_id]

        # Save reference structure (no transformation)
        self._save_split_structure(
            ref_structure,
            ref_dir / f"{ref_pdb_id}_antigen.cif",
            ref_dir / f"{ref_pdb_id}_antibody.cif"
        )

        # Align other members
        results = []

        # Add reference result
        results.append(AlignmentResult(
            epitope_id=ref_epitope_id,
            pdb_id=ref_pdb_id,
            rmsd=0.0,
            num_atoms_aligned=ref_epitope.total_residue_count(),
            rotation_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation_vector=[0, 0, 0],
            is_reference=True
        ))

        rmsds = []

        for member in group.members:
            if member.is_reference:
                continue

            mob_epitope_id = member.epitope_id
            mob_pdb_id = member.pdb_id

            if mob_pdb_id not in structures:
                logger.warning(f"Structure {mob_pdb_id} not found, skipping")
                continue
            if mob_epitope_id not in epitopes:
                logger.warning(f"Epitope {mob_epitope_id} not found, skipping")
                continue

            mob_structure = structures[mob_pdb_id]
            mob_epitope = epitopes[mob_epitope_id]

            try:
                # Perform epitope-based alignment
                rmsd, rotran, num_atoms = self._align_epitope(
                    ref_structure, ref_epitope,
                    mob_structure, mob_epitope
                )

                # Save aligned structure with transformation
                self._save_split_structure_transformed(
                    mob_structure,
                    aligned_dir / f"{mob_pdb_id}_antigen.cif",
                    aligned_dir / f"{mob_pdb_id}_antibody.cif",
                    rotran
                )

                rot, tran = rotran
                results.append(AlignmentResult(
                    epitope_id=mob_epitope_id,
                    pdb_id=mob_pdb_id,
                    rmsd=rmsd,
                    num_atoms_aligned=num_atoms,
                    rotation_matrix=rot.tolist() if hasattr(rot, 'tolist') else list(rot),
                    translation_vector=tran.tolist() if hasattr(tran, 'tolist') else list(tran),
                    is_reference=False
                ))

                rmsds.append(rmsd)
                logger.info(f"Aligned {mob_pdb_id} to {ref_pdb_id}: RMSD={rmsd:.3f}Å")

            except Exception as e:
                logger.error(f"Failed to align {mob_pdb_id}: {e}")
                continue

        avg_rmsd = np.mean(rmsds) if rmsds else 0.0

        output = GroupAlignmentOutput(
            group_id=group.group_id,
            reference_pdb=ref_pdb_id,
            reference_epitope_id=ref_epitope_id,
            members=results,
            avg_rmsd=float(avg_rmsd),
            output_dir=group_dir
        )

        # Save metadata
        self._save_group_metadata(output, group_dir / "group_metadata.json")

        return output

    def _align_epitope(
        self,
        ref_structure: CleanedStructure,
        ref_epitope: EpitopeResidues,
        mob_structure: CleanedStructure,
        mob_epitope: EpitopeResidues
    ) -> Tuple[float, Tuple, int]:
        """
        Align mobile structure to reference based on epitope Cα atoms.

        Args:
            ref_structure: Reference cleaned structure
            ref_epitope: Reference epitope residues
            mob_structure: Mobile cleaned structure
            mob_epitope: Mobile epitope residues

        Returns:
            (rmsd, (rotation, translation), num_atoms_aligned)
        """
        if PYMOL_AVAILABLE:
            return self._align_epitope_pymol(
                ref_structure, ref_epitope,
                mob_structure, mob_epitope
            )
        elif self.fallback_to_biopython:
            return self._align_epitope_biopython(
                ref_structure, ref_epitope,
                mob_structure, mob_epitope
            )
        else:
            raise AlignmentError("PyMOL not available and fallback disabled")

    def _align_epitope_pymol(
        self,
        ref_structure: CleanedStructure,
        ref_epitope: EpitopeResidues,
        mob_structure: CleanedStructure,
        mob_epitope: EpitopeResidues
    ) -> Tuple[float, Tuple, int]:
        """
        Align using PyMOL with epitope-specific selection.
        """
        # Build epitope selection strings for PyMOL
        ref_sel = self._build_pymol_epitope_selection("reference", ref_epitope)
        mob_sel = self._build_pymol_epitope_selection("mobile", mob_epitope)

        # Create temp files
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy source files to temp
            ref_file = tmpdir / "reference.cif"
            mob_file = tmpdir / "mobile.cif"

            import shutil
            shutil.copy(ref_structure.file_path, ref_file)
            shutil.copy(mob_structure.file_path, mob_file)

            # PyMOL script
            align_cmd = "super" if self.use_super else "align"

            script = f'''
import sys
import json
import numpy as np
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-qc'])

cmd.load("{ref_file}", "reference")
cmd.load("{mob_file}", "mobile")

# Align based on epitope Cα atoms
result = cmd.{align_cmd}("{mob_sel} and name CA", "{ref_sel} and name CA")

if isinstance(result, tuple):
    rmsd = result[0]
    n_atoms = result[1]
else:
    rmsd = result
    n_atoms = 0

# Get transformation matrix
# PyMOL stores it internally, we need to extract it
# Using get_object_matrix which returns 4x4 transformation matrix
matrix = cmd.get_object_matrix("mobile")

# Output results
output = {{
    "rmsd": rmsd,
    "n_atoms": n_atoms,
    "matrix": matrix
}}
sys.stderr.write(json.dumps(output) + "\\n")
sys.stderr.flush()
'''

            # Run PyMOL
            result = subprocess.run(
                [PYMOL_PYTHON_PATH, '-c', script],
                capture_output=True,
                text=True,
                timeout=120
            )

            if result.returncode != 0:
                raise AlignmentError(f"PyMOL alignment failed: {result.stderr}")

            # Parse output
            output = None
            for line in result.stderr.strip().split('\n'):
                line = line.strip()
                if line.startswith('{') and 'rmsd' in line:
                    output = json.loads(line)
                    break

            if output is None:
                raise AlignmentError(f"Could not parse PyMOL output: {result.stderr}")

            rmsd = output['rmsd']
            n_atoms = output['n_atoms']

            # Parse 4x4 transformation matrix to rotation + translation
            matrix = output['matrix']
            if matrix and len(matrix) == 16:
                # PyMOL returns column-major 4x4 matrix
                rot = np.array([
                    [matrix[0], matrix[4], matrix[8]],
                    [matrix[1], matrix[5], matrix[9]],
                    [matrix[2], matrix[6], matrix[10]]
                ])
                tran = np.array([matrix[12], matrix[13], matrix[14]])
            else:
                # Fallback: identity
                rot = np.eye(3)
                tran = np.zeros(3)

            return rmsd, (rot, tran), n_atoms

    def _align_epitope_biopython(
        self,
        ref_structure: CleanedStructure,
        ref_epitope: EpitopeResidues,
        mob_structure: CleanedStructure,
        mob_epitope: EpitopeResidues
    ) -> Tuple[float, Tuple, int]:
        """
        Align using Biopython Superimposer with epitope Cα atoms.
        """
        from Bio.PDB import MMCIFParser, Superimposer

        parser = MMCIFParser(QUIET=True)

        ref_bio = parser.get_structure("ref", str(ref_structure.file_path))
        mob_bio = parser.get_structure("mob", str(mob_structure.file_path))

        # Collect epitope Cα atoms
        ref_atoms = []
        mob_atoms = []

        ref_model = ref_bio[0]
        mob_model = mob_bio[0]

        # Get epitope residue indices
        for chain_id, residue_ids in ref_epitope.antigen_chains.items():
            if chain_id not in ref_model:
                continue
            chain = ref_model[chain_id]
            for res_id in residue_ids:
                for residue in chain:
                    if residue.id[1] == res_id and 'CA' in residue:
                        ref_atoms.append(residue['CA'])
                        break

        for chain_id, residue_ids in mob_epitope.antigen_chains.items():
            if chain_id not in mob_model:
                continue
            chain = mob_model[chain_id]
            for res_id in residue_ids:
                for residue in chain:
                    if residue.id[1] == res_id and 'CA' in residue:
                        mob_atoms.append(residue['CA'])
                        break

        # Match by position (assume same order)
        n_atoms = min(len(ref_atoms), len(mob_atoms))
        if n_atoms < 3:
            raise AlignmentError(f"Not enough Cα atoms for alignment: ref={len(ref_atoms)}, mob={len(mob_atoms)}")

        ref_atoms = ref_atoms[:n_atoms]
        mob_atoms = mob_atoms[:n_atoms]

        # Superimpose
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, mob_atoms)

        rmsd = super_imposer.rms
        rot, tran = super_imposer.rotran

        return rmsd, (rot, tran), n_atoms

    def _build_pymol_epitope_selection(
        self,
        object_name: str,
        epitope: EpitopeResidues
    ) -> str:
        """
        Build PyMOL selection string for epitope residues.

        Args:
            object_name: PyMOL object name
            epitope: Epitope residues

        Returns:
            PyMOL selection string like "obj and (chain A and resi 100+101+102)"
        """
        parts = []
        for chain_id, residue_ids in epitope.antigen_chains.items():
            if not residue_ids:
                continue
            resi_str = '+'.join(str(r) for r in residue_ids)
            parts.append(f"(chain {chain_id} and resi {resi_str})")

        if not parts:
            return f"{object_name}"

        return f"{object_name} and ({' or '.join(parts)})"

    def _save_split_structure(
        self,
        structure: CleanedStructure,
        antigen_path: Path,
        antibody_path: Path
    ):
        """
        Save structure as separate antigen and antibody CIF files.

        Args:
            structure: Cleaned structure
            antigen_path: Output path for antigen
            antibody_path: Output path for antibody
        """
        # Get chain IDs
        antigen_chains = [m.original_chain_id for m in structure.chain_mappings
                         if m.chain_type == 'antigen']
        antibody_chains = [m.original_chain_id for m in structure.chain_mappings
                          if m.chain_type in ('antibody_heavy', 'antibody_light', 'antibody')]

        # Save antigen
        self._save_chains_to_cif(
            structure.file_path,
            antigen_path,
            antigen_chains
        )

        # Save antibody
        self._save_chains_to_cif(
            structure.file_path,
            antibody_path,
            antibody_chains
        )

    def _save_split_structure_transformed(
        self,
        structure: CleanedStructure,
        antigen_path: Path,
        antibody_path: Path,
        rotran: Tuple
    ):
        """
        Save structure with transformation applied.

        Args:
            structure: Cleaned structure
            antigen_path: Output path for antigen
            antibody_path: Output path for antibody
            rotran: (rotation_matrix, translation_vector) tuple
        """
        # Get chain IDs
        antigen_chains = [m.original_chain_id for m in structure.chain_mappings
                         if m.chain_type == 'antigen']
        antibody_chains = [m.original_chain_id for m in structure.chain_mappings
                          if m.chain_type in ('antibody_heavy', 'antibody_light', 'antibody')]

        # Save antigen with transformation
        self._save_chains_to_cif_transformed(
            structure.file_path,
            antigen_path,
            antigen_chains,
            rotran
        )

        # Save antibody with transformation
        self._save_chains_to_cif_transformed(
            structure.file_path,
            antibody_path,
            antibody_chains,
            rotran
        )

    def _save_chains_to_cif(
        self,
        source_cif: Path,
        output_path: Path,
        chain_ids: List[str]
    ):
        """
        Save specific chains from CIF file.
        """
        if not chain_ids:
            logger.warning(f"No chains to save for {output_path}")
            return

        chain_set = set(chain_ids)
        doc = gemmi.cif.read_file(str(source_cif))
        block = doc.sole_block()

        table = block.find_mmcif_category("_atom_site")
        if not table:
            doc.write_file(str(output_path))
            return

        tags = list(table.tags)
        idx_auth = tags.index("_atom_site.auth_asym_id") if "_atom_site.auth_asym_id" in tags else -1

        rows_to_remove = []
        for i, row in enumerate(table):
            auth_chain = row[idx_auth] if idx_auth >= 0 else ""
            if auth_chain not in chain_set:
                rows_to_remove.append(i)

        for offset, row_idx in enumerate(rows_to_remove):
            table.remove_row(row_idx - offset)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        doc.write_file(str(output_path))

    def _save_chains_to_cif_transformed(
        self,
        source_cif: Path,
        output_path: Path,
        chain_ids: List[str],
        rotran: Tuple
    ):
        """
        Save specific chains from CIF file with transformation applied.
        """
        if not chain_ids:
            logger.warning(f"No chains to save for {output_path}")
            return

        rot, tran = rotran
        chain_set = set(chain_ids)
        doc = gemmi.cif.read_file(str(source_cif))
        block = doc.sole_block()

        table = block.find_mmcif_category("_atom_site")
        if not table:
            doc.write_file(str(output_path))
            return

        tags = list(table.tags)
        idx_auth = tags.index("_atom_site.auth_asym_id") if "_atom_site.auth_asym_id" in tags else -1
        idx_x = tags.index("_atom_site.Cartn_x") if "_atom_site.Cartn_x" in tags else -1
        idx_y = tags.index("_atom_site.Cartn_y") if "_atom_site.Cartn_y" in tags else -1
        idx_z = tags.index("_atom_site.Cartn_z") if "_atom_site.Cartn_z" in tags else -1

        rows_to_remove = []
        for i, row in enumerate(table):
            auth_chain = row[idx_auth] if idx_auth >= 0 else ""

            if auth_chain not in chain_set:
                rows_to_remove.append(i)
                continue

            # Apply transformation
            if idx_x >= 0 and idx_y >= 0 and idx_z >= 0:
                try:
                    x = float(row[idx_x])
                    y = float(row[idx_y])
                    z = float(row[idx_z])

                    # Apply rotation and translation
                    new_x = rot[0, 0] * x + rot[0, 1] * y + rot[0, 2] * z + tran[0]
                    new_y = rot[1, 0] * x + rot[1, 1] * y + rot[1, 2] * z + tran[1]
                    new_z = rot[2, 0] * x + rot[2, 1] * y + rot[2, 2] * z + tran[2]

                    row[idx_x] = f"{new_x:.3f}"
                    row[idx_y] = f"{new_y:.3f}"
                    row[idx_z] = f"{new_z:.3f}"
                except ValueError:
                    pass

        for offset, row_idx in enumerate(rows_to_remove):
            table.remove_row(row_idx - offset)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        doc.write_file(str(output_path))

    def _save_group_metadata(
        self,
        output: GroupAlignmentOutput,
        output_path: Path
    ):
        """
        Save group alignment metadata to JSON.
        """
        with open(output_path, 'w') as f:
            json.dump(output.to_dict(), f, indent=2)

        logger.info(f"Saved group metadata to {output_path}")


def save_alignment_summary_csv(
    outputs: List[GroupAlignmentOutput],
    output_path: Path
):
    """
    Save alignment summary to CSV.

    Args:
        outputs: List of GroupAlignmentOutput objects
        output_path: Output CSV path
    """
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'group_id', 'reference_pdb', 'member_count', 'avg_rmsd',
            'min_rmsd', 'max_rmsd'
        ])

        for output in outputs:
            rmsds = [m.rmsd for m in output.members if not m.is_reference]

            writer.writerow([
                output.group_id,
                output.reference_pdb,
                len(output.members),
                f"{output.avg_rmsd:.3f}" if rmsds else "N/A",
                f"{min(rmsds):.3f}" if rmsds else "N/A",
                f"{max(rmsds):.3f}" if rmsds else "N/A"
            ])

    logger.info(f"Saved alignment summary to {output_path}")


def generate_alignment_report(outputs: List[GroupAlignmentOutput]) -> str:
    """
    Generate human-readable alignment report.

    Args:
        outputs: List of GroupAlignmentOutput objects

    Returns:
        Formatted report string
    """
    lines = []
    lines.append("=" * 70)
    lines.append("Epitope Alignment Report")
    lines.append("=" * 70)

    total_groups = len(outputs)
    total_structures = sum(len(o.members) for o in outputs)
    all_rmsds = [m.rmsd for o in outputs for m in o.members if not m.is_reference]

    lines.append(f"\nTotal groups aligned: {total_groups}")
    lines.append(f"Total structures: {total_structures}")

    if all_rmsds:
        lines.append(f"\nRMSD statistics:")
        lines.append(f"  Mean: {np.mean(all_rmsds):.3f} Å")
        lines.append(f"  Median: {np.median(all_rmsds):.3f} Å")
        lines.append(f"  Range: {min(all_rmsds):.3f} - {max(all_rmsds):.3f} Å")

    lines.append(f"\nPer-group summary:")
    for output in outputs[:10]:  # Show first 10
        rmsds = [m.rmsd for m in output.members if not m.is_reference]
        rmsd_str = f"avg={output.avg_rmsd:.3f}Å" if rmsds else "reference only"
        lines.append(f"  {output.group_id}: {len(output.members)} members, {rmsd_str}")

    if len(outputs) > 10:
        lines.append(f"  ... and {len(outputs) - 10} more groups")

    lines.append("=" * 70)

    return '\n'.join(lines)
