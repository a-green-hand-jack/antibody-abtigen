#!/usr/bin/env python3
"""
Test script for PyMOLStructureAligner.

Tests:
1. Load groups from test groups.json
2. Load cleaned structures
3. Align structures within groups
4. Verify output files and metadata
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import json
import logging

from antibody_abtigen.epitope_pipeline import (
    GemmiStructureCleaner,
    GeometricEpitopeExtractor,
    PyMOLStructureAligner,
    GroupResult,
    GroupMember,
    save_alignment_summary_csv,
    generate_alignment_report,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    # Paths
    groups_path = project_root / "data/epitope_pipeline/grouping/groups.json"
    structures_path = project_root / "data/cleaned_cif"
    sabdab_path = project_root / "data/meta/sabdab_summary_all.tsv"
    output_dir = project_root / "data/epitope_pipeline/aligned_test"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print("Testing PyMOLStructureAligner")
    print(f"{'='*70}\n")

    # Check if groups.json exists
    if not groups_path.exists():
        print(f"Error: Groups file not found: {groups_path}")
        print("Run 'uv run antibody-abtigen group ...' first to generate groups.json")
        return 1

    # Check if structures exist
    if not structures_path.exists():
        print(f"Error: Structures directory not found: {structures_path}")
        return 1

    # Step 1: Load groups
    print("Step 1: Loading groups...")
    with open(groups_path, 'r') as f:
        groups_data = json.load(f)

    groups_list = groups_data.get('groups', [])
    print(f"  Found {len(groups_list)} groups")

    if not groups_list:
        print("  No groups found. Need at least one group for testing.")
        return 1

    # Take first group for testing
    test_group_data = groups_list[0]
    print(f"  Testing with group: {test_group_data['group_id']}")
    print(f"  Reference: {test_group_data['reference_epitope_id']}")
    print(f"  Members: {test_group_data['member_count']}")

    # Step 2: Initialize components
    print("\nStep 2: Initializing components...")
    cleaner = GemmiStructureCleaner(
        sabdab_summary_path=sabdab_path if sabdab_path.exists() else None
    )
    extractor = GeometricEpitopeExtractor(distance_threshold=5.0)
    aligner = PyMOLStructureAligner(use_super=True)

    # Step 3: Load structures and extract epitopes
    print("\nStep 3: Loading structures and extracting epitopes...")

    pdb_ids = set(m['pdb_id'] for m in test_group_data.get('members', []))
    print(f"  PDB IDs needed: {pdb_ids}")

    # Find CIF files
    cif_files = {f.stem.replace('_cleaned', '').upper(): f
                 for f in structures_path.glob("*_cleaned.cif")}
    print(f"  Found {len(cif_files)} cleaned CIF files")

    structures = {}
    epitopes = {}
    errors = []

    for pdb_id in pdb_ids:
        if pdb_id not in cif_files:
            errors.append(f"CIF file not found for {pdb_id}")
            continue

        cif_path = cif_files[pdb_id]
        try:
            cleaned, _ = cleaner.clean_structure(cif_path, skip_filtering=True)
            if cleaned:
                structures[pdb_id] = cleaned
                print(f"    Loaded {pdb_id}: {len(cleaned.chain_mappings)} chains")

                # Extract epitope
                epitope = extractor.extract_epitope(cleaned)
                epitopes[epitope.epitope_id] = epitope
                print(f"      Epitope: {epitope.total_residue_count()} residues")
        except Exception as e:
            errors.append(f"Failed to process {pdb_id}: {e}")
            print(f"    ERROR: {pdb_id}: {e}")

    print(f"\n  Loaded {len(structures)} structures, {len(epitopes)} epitopes")

    if len(structures) < 2:
        print("  Need at least 2 structures for alignment testing")
        print("  Skipping alignment test")
        return 1

    # Step 4: Convert group data to GroupResult
    print("\nStep 4: Converting group data...")
    members = []
    for m in test_group_data.get('members', []):
        if m['pdb_id'] in structures:
            members.append(GroupMember(
                epitope_id=m['epitope_id'],
                pdb_id=m['pdb_id'],
                is_reference=m['is_reference'],
                antigen_chains=m.get('antigen_chains', []),
                epitope_residues=m.get('epitope_residues', {}),
                total_residues=m.get('total_residues', 0),
                similarity_to_ref=m.get('similarity_to_ref')
            ))

    group = GroupResult(
        group_id=test_group_data['group_id'],
        reference_epitope_id=test_group_data['reference_epitope_id'],
        member_count=len(members),
        avg_similarity=test_group_data['avg_similarity'],
        min_similarity=test_group_data.get('min_similarity', 0),
        max_similarity=test_group_data.get('max_similarity', 1),
        members=members
    )

    print(f"  Group has {len(members)} members (after filtering)")

    # Step 5: Perform alignment
    print("\nStep 5: Performing alignment...")
    try:
        output = aligner.align_group_detailed(group, epitopes, structures, output_dir)

        print(f"\n  Alignment completed!")
        print(f"  Group ID: {output.group_id}")
        print(f"  Reference: {output.reference_pdb}")
        print(f"  Average RMSD: {output.avg_rmsd:.3f} Å")
        print(f"  Members aligned: {len(output.members)}")

        for member in output.members:
            status = "REFERENCE" if member.is_reference else f"RMSD={member.rmsd:.3f}Å"
            print(f"    {member.pdb_id}: {status} ({member.num_atoms_aligned} atoms)")

    except Exception as e:
        print(f"  ERROR: Alignment failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Step 6: Verify output files
    print("\nStep 6: Verifying output files...")
    group_dir = output_dir / output.group_id

    expected_files = [
        group_dir / "group_metadata.json",
        group_dir / "reference" / f"{output.reference_pdb}_antigen.cif",
        group_dir / "reference" / f"{output.reference_pdb}_antibody.cif",
    ]

    for member in output.members:
        if not member.is_reference:
            expected_files.append(group_dir / "aligned" / f"{member.pdb_id}_antigen.cif")
            expected_files.append(group_dir / "aligned" / f"{member.pdb_id}_antibody.cif")

    all_exist = True
    for path in expected_files:
        if path.exists():
            print(f"    ✓ {path.name}")
        else:
            print(f"    ✗ {path.name} MISSING")
            all_exist = False

    # Step 7: Print summary
    print("\n" + "="*70)
    if all_exist:
        print("Aligner test PASSED!")
    else:
        print("Aligner test FAILED - some output files missing")
    print("="*70)
    print(f"Output directory: {output_dir}")

    if errors:
        print(f"\nWarnings ({len(errors)}):")
        for err in errors:
            print(f"  - {err}")

    return 0 if all_exist else 1


if __name__ == "__main__":
    sys.exit(main())
