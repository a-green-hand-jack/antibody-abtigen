"""
Epitope consistency validation module.

This module handles:
1. Extracting epitope residues (antigen residues in contact with antibody)
2. Calculating epitope sequence identity between human and mouse
3. Calculating epitope structure RMSD
4. Validating epitope consistency

An epitope is defined as antigen residues within a distance threshold
of any antibody atom.
"""

import os
from typing import Optional, List, Tuple, Dict, Set
from dataclasses import dataclass

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner


@dataclass
class EpitopeValidationResult:
    """Result of epitope consistency validation."""
    is_consistent: bool
    epitope_identity: float  # Sequence identity percentage (0-100)
    epitope_rmsd: float  # Structure RMSD in Angstroms
    human_epitope_residues: List[Tuple[str, int]]  # (chain_id, residue_num)
    mouse_epitope_residues: List[Tuple[str, int]]
    num_contacts: int  # Number of antibody-antigen contacts
    message: str


def parse_structure_file(file_path: str, structure_id: str = "structure") -> Optional[Structure]:
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

        return parser.get_structure(structure_id, file_path)
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None


def get_epitope_residues(
    structure: Structure,
    antigen_chain_ids: List[str],
    antibody_chain_ids: List[str],
    distance_threshold: float = 5.0
) -> List[Tuple[str, int, str]]:
    """
    Extract epitope residues from a structure.

    Epitope residues are antigen residues that have at least one atom
    within the distance threshold of any antibody atom.

    Args:
        structure: Biopython Structure object
        antigen_chain_ids: List of antigen chain IDs
        antibody_chain_ids: List of antibody chain IDs
        distance_threshold: Distance cutoff in Angstroms (default: 5.0Å)

    Returns:
        List of tuples (chain_id, residue_number, residue_name)
    """
    model = structure[0]

    # Collect all antibody atoms
    antibody_atoms = []
    for chain_id in antibody_chain_ids:
        if chain_id in model:
            for residue in model[chain_id]:
                # Skip water and non-standard residues
                if residue.get_id()[0] != ' ':
                    continue
                for atom in residue:
                    antibody_atoms.append(atom)

    if not antibody_atoms:
        return []

    # Build neighbor search from antibody atoms
    ns = NeighborSearch(antibody_atoms)

    # Find antigen residues in contact with antibody
    epitope_residues: Set[Tuple[str, int, str]] = set()

    for chain_id in antigen_chain_ids:
        if chain_id not in model:
            continue
        for residue in model[chain_id]:
            # Skip water and non-standard residues
            if residue.get_id()[0] != ' ':
                continue

            # Check if any atom in this residue is near an antibody atom
            for atom in residue:
                nearby = ns.search(atom.get_coord(), distance_threshold)
                if nearby:
                    res_num = residue.get_id()[1]
                    res_name = residue.get_resname()
                    epitope_residues.add((chain_id, res_num, res_name))
                    break  # Found contact, move to next residue

    return sorted(list(epitope_residues), key=lambda x: (x[0], x[1]))


def extract_epitope_sequence(
    structure: Structure,
    epitope_residues: List[Tuple[str, int, str]]
) -> str:
    """
    Extract the amino acid sequence of epitope residues.

    Args:
        structure: Biopython Structure object
        epitope_residues: List of (chain_id, residue_num, residue_name) tuples

    Returns:
        One-letter amino acid sequence string
    """
    sequence = []
    for chain_id, res_num, res_name in epitope_residues:
        try:
            aa = seq1(res_name)
            if aa != 'X':  # Skip unknown residues
                sequence.append(aa)
        except Exception:
            pass

    return ''.join(sequence)


def calculate_epitope_sequence_identity(
    human_epitope_seq: str,
    mouse_epitope_seq: str
) -> float:
    """
    Calculate sequence identity between human and mouse epitope sequences.

    Args:
        human_epitope_seq: Human epitope sequence
        mouse_epitope_seq: Mouse epitope sequence

    Returns:
        Sequence identity as percentage (0-100)
    """
    if not human_epitope_seq or not mouse_epitope_seq:
        return 0.0

    # Perform global alignment
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(human_epitope_seq, mouse_epitope_seq))
    if not alignments:
        return 0.0

    alignment = alignments[0]

    # Count matches
    matches = 0
    aligned_human = alignment[0]
    aligned_mouse = alignment[1]

    for i in range(len(aligned_human)):
        if aligned_human[i] != '-' and aligned_mouse[i] != '-':
            if aligned_human[i] == aligned_mouse[i]:
                matches += 1

    # Calculate identity as matches / length of shorter sequence
    shorter_len = min(len(human_epitope_seq), len(mouse_epitope_seq))
    if shorter_len == 0:
        return 0.0

    return (matches / shorter_len) * 100.0


def extract_epitope_ca_coords(
    structure: Structure,
    epitope_residues: List[Tuple[str, int, str]]
) -> np.ndarray:
    """
    Extract CA atom coordinates for epitope residues.

    Args:
        structure: Biopython Structure object
        epitope_residues: List of (chain_id, residue_num, residue_name) tuples

    Returns:
        Numpy array of CA coordinates (N x 3)
    """
    model = structure[0]
    coords = []

    for chain_id, res_num, _ in epitope_residues:
        if chain_id not in model:
            continue
        chain = model[chain_id]

        for residue in chain:
            if residue.get_id()[1] == res_num and residue.get_id()[0] == ' ':
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
                break

    return np.array(coords)


def calculate_epitope_rmsd(
    human_structure: Structure,
    mouse_structure: Structure,
    human_epitope_residues: List[Tuple[str, int, str]],
    mouse_epitope_residues: List[Tuple[str, int, str]]
) -> float:
    """
    Calculate RMSD between human and mouse epitope regions.

    Note: This assumes the structures are already aligned.

    Args:
        human_structure: Human antigen structure (aligned)
        mouse_structure: Mouse antigen structure (aligned to human)
        human_epitope_residues: Human epitope residue list
        mouse_epitope_residues: Mouse epitope residue list

    Returns:
        RMSD in Angstroms
    """
    human_coords = extract_epitope_ca_coords(human_structure, human_epitope_residues)
    mouse_coords = extract_epitope_ca_coords(mouse_structure, mouse_epitope_residues)

    if len(human_coords) == 0 or len(mouse_coords) == 0:
        return float('inf')

    # Use the shorter list for comparison
    min_len = min(len(human_coords), len(mouse_coords))
    human_coords = human_coords[:min_len]
    mouse_coords = mouse_coords[:min_len]

    # Calculate RMSD
    diff = human_coords - mouse_coords
    rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

    return rmsd


def validate_epitope_consistency(
    human_complex_file: str,
    mouse_antigen_file: str,
    human_antigen_chain: str,
    mouse_antigen_chain: str,
    antibody_chains: List[str],
    distance_threshold: float = 5.0,
    rmsd_threshold: float = 1.5,
    identity_threshold: float = 80.0
) -> EpitopeValidationResult:
    """
    Validate epitope consistency between human and mouse antigens.

    This function:
    1. Extracts epitope residues from the human complex (antigen + antibody)
    2. Maps epitope residues to mouse antigen by sequence alignment
    3. Calculates sequence identity and structure RMSD
    4. Determines if the epitope is consistent

    Args:
        human_complex_file: Path to human complex structure (with antibody)
        mouse_antigen_file: Path to aligned mouse antigen structure
        human_antigen_chain: Human antigen chain ID
        mouse_antigen_chain: Mouse antigen chain ID
        antibody_chains: List of antibody chain IDs
        distance_threshold: Contact distance threshold (default: 5.0Å)
        rmsd_threshold: Maximum allowed epitope RMSD (default: 1.5Å)
        identity_threshold: Minimum required epitope identity (default: 80%)

    Returns:
        EpitopeValidationResult with validation status and metrics
    """
    # Parse structures
    human_structure = parse_structure_file(human_complex_file, "human")
    mouse_structure = parse_structure_file(mouse_antigen_file, "mouse")

    if human_structure is None:
        return EpitopeValidationResult(
            is_consistent=False,
            epitope_identity=0.0,
            epitope_rmsd=float('inf'),
            human_epitope_residues=[],
            mouse_epitope_residues=[],
            num_contacts=0,
            message=f"Could not parse human structure: {human_complex_file}"
        )

    if mouse_structure is None:
        return EpitopeValidationResult(
            is_consistent=False,
            epitope_identity=0.0,
            epitope_rmsd=float('inf'),
            human_epitope_residues=[],
            mouse_epitope_residues=[],
            num_contacts=0,
            message=f"Could not parse mouse structure: {mouse_antigen_file}"
        )

    # Extract human epitope residues
    human_epitope = get_epitope_residues(
        human_structure,
        [human_antigen_chain],
        antibody_chains,
        distance_threshold
    )

    if not human_epitope:
        return EpitopeValidationResult(
            is_consistent=False,
            epitope_identity=0.0,
            epitope_rmsd=float('inf'),
            human_epitope_residues=[],
            mouse_epitope_residues=[],
            num_contacts=0,
            message="No epitope residues found in human structure"
        )

    # Extract sequences
    human_epitope_seq = extract_epitope_sequence(human_structure, human_epitope)

    # For mouse, we need to find corresponding residues
    # First get all residues in mouse antigen chain
    mouse_model = mouse_structure[0]
    mouse_residues = []
    if mouse_antigen_chain in mouse_model:
        for residue in mouse_model[mouse_antigen_chain]:
            if residue.get_id()[0] == ' ':  # Standard residue
                mouse_residues.append((
                    mouse_antigen_chain,
                    residue.get_id()[1],
                    residue.get_resname()
                ))

    # Extract mouse epitope based on aligned positions
    # Use sequence alignment to find corresponding residues
    human_full_seq = ""
    human_res_map: Dict[int, int] = {}  # seq_pos -> res_num
    pos = 0
    for chain_id, res_num, res_name in sorted(
        [(human_antigen_chain, r.get_id()[1], r.get_resname())
         for r in human_structure[0][human_antigen_chain] if r.get_id()[0] == ' '],
        key=lambda x: x[1]
    ):
        aa = seq1(res_name)
        if aa != 'X':
            human_full_seq += aa
            human_res_map[pos] = res_num
            pos += 1

    mouse_full_seq = ""
    mouse_res_map: Dict[int, int] = {}  # seq_pos -> res_num
    pos = 0
    for chain_id, res_num, res_name in mouse_residues:
        aa = seq1(res_name)
        if aa != 'X':
            mouse_full_seq += aa
            mouse_res_map[pos] = res_num
            pos += 1

    # Align sequences to find corresponding positions
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(human_full_seq, mouse_full_seq))
    if not alignments:
        return EpitopeValidationResult(
            is_consistent=False,
            epitope_identity=0.0,
            epitope_rmsd=float('inf'),
            human_epitope_residues=[(c, n) for c, n, _ in human_epitope],
            mouse_epitope_residues=[],
            num_contacts=len(human_epitope),
            message="Could not align human and mouse sequences"
        )

    alignment = alignments[0]

    # Build residue number mapping
    human_to_mouse: Dict[int, int] = {}  # human_res_num -> mouse_res_num
    h_idx = 0
    m_idx = 0

    aligned_human = alignment[0]
    aligned_mouse = alignment[1]

    for i in range(len(aligned_human)):
        h_char = aligned_human[i]
        m_char = aligned_mouse[i]

        if h_char != '-' and m_char != '-':
            if h_idx in human_res_map and m_idx in mouse_res_map:
                human_to_mouse[human_res_map[h_idx]] = mouse_res_map[m_idx]

        if h_char != '-':
            h_idx += 1
        if m_char != '-':
            m_idx += 1

    # Map human epitope to mouse
    mouse_epitope = []
    for chain_id, res_num, res_name in human_epitope:
        if res_num in human_to_mouse:
            mouse_res_num = human_to_mouse[res_num]
            # Find mouse residue name
            for mc, mn, mr in mouse_residues:
                if mn == mouse_res_num:
                    mouse_epitope.append((mouse_antigen_chain, mouse_res_num, mr))
                    break

    mouse_epitope_seq = extract_epitope_sequence(mouse_structure, mouse_epitope)

    # Calculate sequence identity
    epitope_identity = calculate_epitope_sequence_identity(
        human_epitope_seq, mouse_epitope_seq
    )

    # Calculate RMSD
    epitope_rmsd = calculate_epitope_rmsd(
        human_structure, mouse_structure,
        human_epitope, mouse_epitope
    )

    # Determine consistency
    is_consistent = (
        epitope_identity >= identity_threshold and
        epitope_rmsd <= rmsd_threshold
    )

    message = f"Epitope identity: {epitope_identity:.1f}% (threshold: {identity_threshold}%), "
    message += f"RMSD: {epitope_rmsd:.2f}Å (threshold: {rmsd_threshold}Å)"

    return EpitopeValidationResult(
        is_consistent=is_consistent,
        epitope_identity=epitope_identity,
        epitope_rmsd=epitope_rmsd,
        human_epitope_residues=[(c, n) for c, n, _ in human_epitope],
        mouse_epitope_residues=[(c, n) for c, n, _ in mouse_epitope],
        num_contacts=len(human_epitope),
        message=message
    )


if __name__ == "__main__":
    # Test the module
    print("Epitope validation module loaded successfully")
