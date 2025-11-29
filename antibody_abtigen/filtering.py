"""
Interaction-based filtering utilities.

Checks whether antibody and antigen chains make contact in saved PDBs.
"""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass
from typing import List, Tuple

from Bio.PDB import NeighborSearch, PDBParser


@dataclass
class FilterResult:
    dp_id: str
    contact_count: int
    passed: bool
    reason: str


def _load_atoms(pdb_path: str):
    """Load atoms from a PDB file; returns a flat list of Atom objects."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    return [atom for atom in structure.get_atoms()]


def _count_contacts(
    antibody_atoms,
    antigen_atoms,
    distance_threshold: float,
    min_contacts: int,
) -> Tuple[int, bool]:
    """Return total contacts and pass/fail."""
    ns = NeighborSearch(list(antibody_atoms))
    contact_count = 0

    for ag_atom in antigen_atoms:
        neighbors = ns.search(ag_atom.get_coord(), distance_threshold)
        contact_count += len(neighbors)
        if contact_count >= min_contacts:
            return contact_count, True

    return contact_count, contact_count >= min_contacts


def evaluate_data_point(
    dp_dir: str,
    distance_threshold: float = 5.0,
    min_contacts: int = 10,
) -> FilterResult:
    """Check if antibody and antigen chains in a DP folder are in contact."""
    metadata_path = os.path.join(dp_dir, "metadata.json")
    if not os.path.exists(metadata_path):
        return FilterResult(
            dp_id=os.path.basename(dp_dir),
            contact_count=0,
            passed=False,
            reason="metadata.json not found",
        )

    try:
        with open(metadata_path, "r") as f:
            metadata = json.load(f)
    except Exception as e:
        return FilterResult(
            dp_id=os.path.basename(dp_dir),
            contact_count=0,
            passed=False,
            reason=f"metadata load error: {e}",
        )

    files = metadata.get("files", {})
    antibody_pdb = files.get("antibody_pdb")
    human_ag_pdb = files.get("human_ag_pdb")

    if not antibody_pdb or not human_ag_pdb:
        return FilterResult(
            dp_id=metadata.get("id", os.path.basename(dp_dir)),
            contact_count=0,
            passed=False,
            reason="missing antibody or antigen file",
        )

    antibody_path = os.path.join(dp_dir, antibody_pdb)
    antigen_path = os.path.join(dp_dir, human_ag_pdb)

    if not os.path.exists(antibody_path) or not os.path.exists(antigen_path):
        return FilterResult(
            dp_id=metadata.get("id", os.path.basename(dp_dir)),
            contact_count=0,
            passed=False,
            reason="PDB files not found",
        )

    try:
        antibody_atoms = _load_atoms(antibody_path)
        antigen_atoms = _load_atoms(antigen_path)
        contact_count, passed = _count_contacts(
            antibody_atoms,
            antigen_atoms,
            distance_threshold,
            min_contacts,
        )
        return FilterResult(
            dp_id=metadata.get("id", os.path.basename(dp_dir)),
            contact_count=contact_count,
            passed=passed,
            reason="" if passed else "no contact",
        )
    except Exception as e:
        return FilterResult(
            dp_id=metadata.get("id", os.path.basename(dp_dir)),
            contact_count=0,
            passed=False,
            reason=f"contact check error: {e}",
        )


def filter_dataset(
    input_dir: str,
    output_dir: str,
    distance_threshold: float = 5.0,
    min_contacts: int = 10,
    dry_run: bool = False,
) -> List[FilterResult]:
    """
    Filter DP_* folders by antibody-antigen contacts.

    Copies passing folders to output_dir unless dry_run is True.
    """
    os.makedirs(output_dir, exist_ok=True)
    results: List[FilterResult] = []

    for name in sorted(os.listdir(input_dir)):
        dp_dir = os.path.join(input_dir, name)
        if not os.path.isdir(dp_dir) or not name.startswith("DP_"):
            continue

        result = evaluate_data_point(
            dp_dir,
            distance_threshold=distance_threshold,
            min_contacts=min_contacts,
        )
        results.append(result)

        if result.passed and not dry_run:
            target = os.path.join(output_dir, name)
            shutil.copytree(dp_dir, target, dirs_exist_ok=True)

    # Write summary CSV
    summary_path = os.path.join(output_dir, "filter_summary.csv")
    try:
        with open(summary_path, "w") as f:
            f.write("id,contact_count,passed,reason\n")
            for r in results:
                f.write(f"{r.dp_id},{r.contact_count},{r.passed},{r.reason}\n")
    except Exception:
        # If writing fails, we still return results
        pass

    return results
