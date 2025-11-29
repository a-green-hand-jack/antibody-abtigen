"""
PDB structure retrieval and processing module.

This module handles:
1. Downloading PDB/mmCIF files
2. Parsing and extracting chains
3. Structure alignment using PyMOL
4. Saving processed structures
"""

import os
import requests
import tempfile
from typing import Optional, List, Dict, Tuple
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# Try to import PyMOL
PYMOL_AVAILABLE = False
PYMOL_APP_PYTHON = "/Applications/PyMOL.app/Contents/bin/python"

try:
    import pymol
    from pymol import cmd
    # Test if it actually works
    pymol.finish_launching(['pymol', '-qc'])
    PYMOL_AVAILABLE = True
except (ImportError, Exception) as e:
    # Check if PyMOL.app is available as fallback
    if os.path.exists(PYMOL_APP_PYTHON):
        PYMOL_AVAILABLE = True
        print("Using PyMOL.app for alignment.")
    else:
        print("Warning: PyMOL not available. Using Biopython for alignment.")


PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
PDB_DOWNLOAD_URL_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"


class ChainSelect(Select):
    """Select specific chains from a structure."""

    def __init__(self, chain_ids: List[str], remove_water: bool = True, remove_hetatm: bool = False):
        self.chain_ids = chain_ids
        self.remove_water = remove_water
        self.remove_hetatm = remove_hetatm

    def accept_chain(self, chain):
        return 1 if chain.get_id() in self.chain_ids else 0

    def accept_residue(self, residue):
        hetflag = residue.get_id()[0]

        # Water molecules
        if self.remove_water and hetflag == 'W':
            return 0

        # HETATM records (non-standard residues)
        if self.remove_hetatm and hetflag.startswith('H_'):
            return 0

        return 1


def download_structure(pdb_id: str, output_dir: str, prefer_cif: bool = True) -> Tuple[Optional[str], Optional[str]]:
    """
    Download a PDB structure file.

    Args:
        pdb_id: PDB ID
        output_dir: Directory to save the file
        prefer_cif: If True, try to download CIF first

    Returns:
        Tuple of (cif_path, pdb_path) - one may be None
    """
    os.makedirs(output_dir, exist_ok=True)

    cif_path = None
    pdb_path = None

    # Download CIF
    cif_file = os.path.join(output_dir, f"{pdb_id.upper()}.cif")
    if not os.path.exists(cif_file):
        try:
            url = PDB_DOWNLOAD_URL.format(pdb_id=pdb_id.upper())
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            with open(cif_file, 'wb') as f:
                f.write(response.content)
            cif_path = cif_file
        except requests.RequestException as e:
            print(f"Could not download CIF for {pdb_id}: {e}")
    else:
        cif_path = cif_file

    # Download PDB format as well
    pdb_file = os.path.join(output_dir, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(pdb_file):
        try:
            url = PDB_DOWNLOAD_URL_PDB.format(pdb_id=pdb_id.upper())
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            with open(pdb_file, 'wb') as f:
                f.write(response.content)
            pdb_path = pdb_file
        except requests.RequestException as e:
            print(f"Could not download PDB for {pdb_id}: {e}")
    else:
        pdb_path = pdb_file

    return cif_path, pdb_path


def parse_structure(file_path: str, structure_id: str = "structure") -> Optional[Structure]:
    """
    Parse a structure file (CIF or PDB).

    Args:
        file_path: Path to structure file
        structure_id: ID to assign to the structure

    Returns:
        Biopython Structure object or None
    """
    if not os.path.exists(file_path):
        return None

    try:
        if file_path.endswith('.cif'):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure(structure_id, file_path)
        return structure

    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None


def extract_chains(
    structure: Structure,
    chain_ids: List[str],
    remove_water: bool = True,
    remove_hetatm: bool = False
) -> Optional[Structure]:
    """
    Extract specific chains from a structure.

    Args:
        structure: Biopython Structure object
        chain_ids: List of chain IDs to extract
        remove_water: Remove water molecules
        remove_hetatm: Remove HETATM records

    Returns:
        New Structure with only specified chains
    """
    from Bio.PDB import StructureBuilder

    builder = StructureBuilder.StructureBuilder()
    builder.init_structure(structure.id + "_extracted")
    builder.init_model(0)

    model = structure[0]

    for chain_id in chain_ids:
        if chain_id in model:
            chain = model[chain_id]
            builder.init_chain(chain_id)

            for residue in chain:
                hetflag = residue.get_id()[0]

                # Skip water if requested
                if remove_water and hetflag == 'W':
                    continue

                # Skip HETATM if requested
                if remove_hetatm and hetflag.startswith('H_'):
                    continue

                builder.init_residue(
                    residue.get_resname(),
                    residue.get_id()[0],
                    residue.get_id()[1],
                    residue.get_id()[2]
                )

                for atom in residue:
                    builder.init_atom(
                        atom.get_name(),
                        atom.get_coord(),
                        atom.get_bfactor(),
                        atom.get_occupancy(),
                        atom.get_altloc(),
                        atom.get_fullname(),
                        atom.get_serial_number(),
                        atom.element
                    )

    return builder.get_structure()


def save_structure(structure: Structure, output_path: str, chain_ids: Optional[List[str]] = None):
    """
    Save a structure to PDB format.

    Args:
        structure: Biopython Structure object
        output_path: Output file path
        chain_ids: Optional list of chains to save (None = all chains)
    """
    io = PDBIO()
    io.set_structure(structure)

    if chain_ids:
        io.save(output_path, ChainSelect(chain_ids))
    else:
        io.save(output_path)


def save_structure_cif(structure: Structure, output_path: str, chain_ids: Optional[List[str]] = None):
    """
    Save a structure to mmCIF format.

    Args:
        structure: Biopython Structure object
        output_path: Output file path
        chain_ids: Optional list of chains to save
    """
    from Bio.PDB.mmcifio import MMCIFIO

    io = MMCIFIO()
    io.set_structure(structure)

    if chain_ids:
        io.save(output_path, ChainSelect(chain_ids))
    else:
        io.save(output_path)


def align_structures_biopython(
    mobile_structure: Structure,
    reference_structure: Structure,
    mobile_chain: str,
    reference_chain: str
) -> Tuple[float, Structure]:
    """
    Align mobile structure to reference using Biopython Superimposer.

    Uses sequence alignment to find corresponding residues before structural alignment.

    Args:
        mobile_structure: Structure to move
        reference_structure: Reference structure (stays fixed)
        mobile_chain: Chain ID in mobile structure
        reference_chain: Chain ID in reference structure

    Returns:
        Tuple of (RMSD, aligned structure)
    """
    from Bio.Align import PairwiseAligner
    from Bio.SeqUtils import seq1

    ref_chain_obj = reference_structure[0][reference_chain]
    mob_chain_obj = mobile_structure[0][mobile_chain]

    # Get residues (exclude water and heteroatoms)
    ref_residues = [r for r in ref_chain_obj if r.get_id()[0] == ' ']
    mob_residues = [r for r in mob_chain_obj if r.get_id()[0] == ' ']

    # Extract sequences
    ref_seq = ''.join([seq1(r.get_resname()) for r in ref_residues])
    mob_seq = ''.join([seq1(r.get_resname()) for r in mob_residues])

    # Perform sequence alignment
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(ref_seq, mob_seq))
    if not alignments:
        raise ValueError("Could not align sequences")

    alignment = alignments[0]

    # Get aligned positions
    ref_atoms = []
    mob_atoms = []

    # Parse alignment to get corresponding residue indices
    ref_idx = 0
    mob_idx = 0

    aligned_ref = alignment[0]
    aligned_mob = alignment[1]

    for i in range(len(aligned_ref)):
        ref_char = aligned_ref[i]
        mob_char = aligned_mob[i]

        if ref_char != '-' and mob_char != '-':
            # Both have residues at this position
            if ref_idx < len(ref_residues) and mob_idx < len(mob_residues):
                ref_res = ref_residues[ref_idx]
                mob_res = mob_residues[mob_idx]

                if 'CA' in ref_res and 'CA' in mob_res:
                    ref_atoms.append(ref_res['CA'])
                    mob_atoms.append(mob_res['CA'])

        if ref_char != '-':
            ref_idx += 1
        if mob_char != '-':
            mob_idx += 1

    if len(ref_atoms) < 3:
        raise ValueError(f"Not enough CA atoms for alignment (found {len(ref_atoms)})")

    # Perform superimposition
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, mob_atoms)

    # Apply transformation to all atoms in mobile structure
    super_imposer.apply(mobile_structure.get_atoms())

    rmsd = super_imposer.rms

    return rmsd, mobile_structure


def align_structures_pymol(
    mobile_file: str,
    reference_file: str,
    output_file: str,
    mobile_chain: Optional[str] = None,
    reference_chain: Optional[str] = None
) -> float:
    """
    Align mobile structure to reference using PyMOL.

    Args:
        mobile_file: Path to mobile structure file
        reference_file: Path to reference structure file
        output_file: Path to save aligned structure
        mobile_chain: Optional chain to align in mobile
        reference_chain: Optional chain to align in reference

    Returns:
        RMSD after alignment
    """
    if not PYMOL_AVAILABLE:
        raise ImportError("PyMOL is not available")

    # Try using PyMOL.app first (more reliable on macOS)
    if os.path.exists(PYMOL_APP_PYTHON):
        return _align_via_pymol_app(
            mobile_file, reference_file, output_file,
            mobile_chain, reference_chain
        )

    # Fallback to direct PyMOL module
    return _align_via_pymol_module(
        mobile_file, reference_file, output_file,
        mobile_chain, reference_chain
    )


def _align_via_pymol_app(
    mobile_file: str,
    reference_file: str,
    output_file: str,
    mobile_chain: Optional[str] = None,
    reference_chain: Optional[str] = None
) -> float:
    """
    Align structures using PyMOL.app's Python interpreter.
    """
    import subprocess
    import json

    # Build selection strings
    ref_sel = "reference"
    mob_sel = "mobile"
    if reference_chain:
        ref_sel = f"reference and chain {reference_chain}"
    if mobile_chain:
        mob_sel = f"mobile and chain {mobile_chain}"

    # Create a PyMOL script
    script = f'''
import sys
import pymol
from pymol import cmd
import json

pymol.finish_launching(['pymol', '-qc'])

cmd.load("{reference_file}", "reference")
cmd.load("{mobile_file}", "mobile")

result = cmd.align("{mob_sel}", "{ref_sel}")
rmsd = result[0] if isinstance(result, tuple) else result

cmd.save("{output_file}", "mobile")

# Output to stderr since PyMOL may capture stdout
sys.stderr.write(json.dumps({{"rmsd": rmsd}}) + "\\n")
sys.stderr.flush()
'''

    # Run via PyMOL.app's Python
    result = subprocess.run(
        [PYMOL_APP_PYTHON, '-c', script],
        capture_output=True,
        text=True,
        timeout=120
    )

    if result.returncode != 0:
        raise RuntimeError(f"PyMOL alignment failed: {result.stderr}")

    # Parse RMSD from stderr output (PyMOL captures stdout)
    output = result.stderr.strip()
    for line in output.split('\n'):
        line = line.strip()
        if line.startswith('{') and 'rmsd' in line:
            data = json.loads(line)
            return data['rmsd']

    raise RuntimeError(f"Could not parse RMSD from PyMOL output. stdout: {result.stdout}, stderr: {result.stderr}")


def _align_via_pymol_module(
    mobile_file: str,
    reference_file: str,
    output_file: str,
    mobile_chain: Optional[str] = None,
    reference_chain: Optional[str] = None
) -> float:
    """
    Align structures using imported PyMOL module.
    """
    # Initialize PyMOL in quiet mode
    pymol.finish_launching(['pymol', '-qc'])

    try:
        # Load structures
        cmd.load(reference_file, "reference")
        cmd.load(mobile_file, "mobile")

        # Build selection strings
        ref_sel = "reference"
        mob_sel = "mobile"

        if reference_chain:
            ref_sel = f"reference and chain {reference_chain}"
        if mobile_chain:
            mob_sel = f"mobile and chain {mobile_chain}"

        # Perform alignment
        result = cmd.align(mob_sel, ref_sel)
        rmsd = result[0] if isinstance(result, tuple) else result

        # Save aligned mobile structure
        cmd.save(output_file, "mobile")

        return rmsd

    finally:
        # Clean up
        cmd.delete("all")


def get_structure_resolution(file_path: str) -> Optional[float]:
    """
    Get resolution from a structure file.

    Args:
        file_path: Path to structure file

    Returns:
        Resolution in Angstroms or None
    """
    try:
        if file_path.endswith('.cif'):
            cif_dict = MMCIF2Dict(file_path)
            resolution = cif_dict.get('_refine.ls_d_res_high', [None])[0]
            if resolution:
                return float(resolution)
        else:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('REMARK   2 RESOLUTION.'):
                        parts = line.split()
                        for i, part in enumerate(parts):
                            if part == 'RESOLUTION.':
                                if i + 1 < len(parts):
                                    try:
                                        return float(parts[i + 1])
                                    except ValueError:
                                        pass
        return None
    except Exception:
        return None


def check_missing_residues(structure: Structure, chain_id: str) -> List[int]:
    """
    Check for missing residues in a chain.

    Args:
        structure: Biopython Structure
        chain_id: Chain to check

    Returns:
        List of missing residue numbers
    """
    chain = structure[0][chain_id]
    residue_nums = [r.get_id()[1] for r in chain if r.get_id()[0] == ' ']

    if not residue_nums:
        return []

    missing = []
    for i in range(min(residue_nums), max(residue_nums) + 1):
        if i not in residue_nums:
            missing.append(i)

    return missing


if __name__ == "__main__":
    # Test the module
    print("Testing structure download...")
    test_dir = tempfile.mkdtemp()

    cif_path, pdb_path = download_structure("5FUO", test_dir)
    print(f"Downloaded: CIF={cif_path}, PDB={pdb_path}")

    if cif_path:
        print("\nTesting structure parsing...")
        structure = parse_structure(cif_path, "5FUO")
        if structure:
            print(f"Parsed structure with {len(list(structure.get_chains()))} chains")

            print("\nChains in structure:")
            for chain in structure[0]:
                print(f"  Chain {chain.id}: {len(list(chain.get_residues()))} residues")
