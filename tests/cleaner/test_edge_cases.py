"""
Test cleaner behavior on edge cases:
1. No antigen (9ia4)
2. Peptide antigen (9ky2)
3. Protein antigen (9jy2) - control
"""

from pathlib import Path
from antibody_abtigen.epitope_pipeline import GemmiStructureCleaner, GeometricEpitopeExtractor

def test_structure(pdb_id, description):
    print(f"\n{'=' * 70}")
    print(f"Testing: {pdb_id} - {description}")
    print('=' * 70)

    raw_cif = Path(f"data/raw_cif/{pdb_id.lower()}.cif")
    if not raw_cif.exists():
        print(f"  ✗ CIF file not found: {raw_cif}")
        return

    output_dir = Path("tests/cleaner/output_edge_cases")
    output_dir.mkdir(parents=True, exist_ok=True)

    sabdab_summary = Path("data/meta/sabdab_summary_all.tsv")
    cleaner = GemmiStructureCleaner(
        output_dir=output_dir,
        sabdab_summary_path=sabdab_summary
    )
    extractor = GeometricEpitopeExtractor(distance_threshold=5.0)

    try:
        # Clean structure
        cleaned = cleaner.clean_structure(raw_cif)

        print(f"\n✓ Cleaned successfully")
        print(f"  PDB ID: {cleaned.pdb_id}")
        print(f"  Total chains: {len(cleaned.chain_mappings)}")

        # Show chain details
        antigen_chains = cleaned.get_antigen_chains()
        antibody_chains = cleaned.get_antibody_chains()

        print(f"\n  Antigen chains: {len(antigen_chains)}")
        for m in antigen_chains:
            print(f"    {m.original_chain_id} ({m.chain_type}): {len(m.sequence)} residues")
            print(f"      Sequence: {m.sequence[:50]}..." if len(m.sequence) > 50 else f"      Sequence: {m.sequence}")

        print(f"\n  Antibody chains: {len(antibody_chains)}")
        for m in antibody_chains:
            print(f"    {m.original_chain_id} ({m.chain_type}): {len(m.sequence)} residues")

        # Try epitope extraction
        if antigen_chains and antibody_chains:
            print(f"\n  Attempting epitope extraction...")
            epitope = extractor.extract_epitope(cleaned)
            print(f"  ✓ Epitope extracted: {epitope.total_residue_count()} residues")
            print(f"    Antigen chains with epitope: {list(epitope.antigen_chains.keys())}")
            for chain_id, residues in epitope.antigen_chains.items():
                print(f"      {chain_id}: {len(residues)} epitope residues")
        else:
            print(f"\n  ⚠ Cannot extract epitope:")
            if not antigen_chains:
                print(f"    - No antigen chains found")
            if not antibody_chains:
                print(f"    - No antibody chains found")

    except Exception as e:
        print(f"\n✗ Failed: {e}")
        import traceback
        traceback.print_exc()

def main():
    print("=" * 70)
    print("Edge Case Testing for Structure Cleaner")
    print("=" * 70)

    # Test cases
    test_structure("9jy2", "Normal protein antigen (control)")
    test_structure("9ky2", "Protein + Peptide antigen")
    test_structure("9ia4", "No antigen (NA)")

if __name__ == "__main__":
    main()
