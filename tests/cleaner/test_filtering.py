"""
Test antigen type filtering functionality.

This script validates the filtering logic on different structure types:
1. Normal protein antigen (should accept)
2. Protein + peptide antigen (should accept - has protein)
3. Peptide-only antigen (should reject)
4. No antigen (should reject)
"""

from pathlib import Path
from antibody_abtigen.epitope_pipeline import (
    GemmiStructureCleaner,
    FilterReason,
    save_filter_log,
    generate_filter_summary,
)

def main():
    print("=" * 70)
    print("Testing Antigen Type Filtering")
    print("=" * 70)

    raw_cif_dir = Path("data/raw_cif")
    output_dir = Path("tests/cleaner/output_filtered")
    output_dir.mkdir(parents=True, exist_ok=True)

    sabdab_summary = Path("data/meta/sabdab_summary_all.tsv")

    # Initialize cleaner
    cleaner = GemmiStructureCleaner(
        output_dir=output_dir,
        sabdab_summary_path=sabdab_summary
    )

    # Test cases
    test_cases = [
        ("9jy2", "Normal protein antigen"),
        ("9ky2", "Protein + peptide antigen"),
        ("9ia3", "Peptide-only antigen"),
        ("9ia4", "No antigen (NA)"),
    ]

    filter_results = []

    for pdb_id, description in test_cases:
        print(f"\n{'=' * 70}")
        print(f"Test: {pdb_id} - {description}")
        print('=' * 70)

        cif_path = raw_cif_dir / f"{pdb_id.lower()}.cif"
        if not cif_path.exists():
            print(f"  ✗ CIF file not found: {cif_path}")
            continue

        try:
            cleaned, filter_result = cleaner.clean_structure(cif_path)
            filter_results.append(filter_result)

            print(f"\n  PDB ID: {filter_result.pdb_id}")
            print(f"  Accepted: {'Yes' if filter_result.accepted else 'No'}")
            print(f"  Reason: {filter_result.reason.value}")
            print(f"  Antigen chains: {filter_result.antigen_chains}")
            print(f"  Antigen type: {filter_result.antigen_type}")
            print(f"  Num antigen chains: {filter_result.num_antigen_chains}")
            print(f"  Total antigen residues: {filter_result.total_antigen_residues}")
            print(f"  Num antibody chains: {filter_result.num_antibody_chains}")
            print(f"  Details: {filter_result.details}")

            if cleaned:
                print(f"\n  ✓ Structure cleaned and saved:")
                print(f"    {cleaned.file_path}")
            else:
                print(f"\n  ⚠ Structure filtered out (not saved)")

        except Exception as e:
            print(f"\n  ✗ Error: {e}")
            import traceback
            traceback.print_exc()

    # Test on larger batch
    print("\n" + "=" * 70)
    print("Batch Testing (First 50 Structures)")
    print("=" * 70)

    cif_files = sorted(raw_cif_dir.glob("*.cif"))[:50]
    batch_results = []

    for cif_path in cif_files:
        try:
            cleaned, filter_result = cleaner.clean_structure(cif_path)
            batch_results.append(filter_result)
        except Exception as e:
            print(f"Warning: Failed to process {cif_path.name}: {e}")

    # Save filter log
    log_path = Path("tests/cleaner/filtering_log.csv")
    save_filter_log(batch_results, log_path)
    print(f"\n✓ Filter log saved to: {log_path}")

    # Generate and print summary
    summary = generate_filter_summary(batch_results)
    print("\n" + summary)

    # Show examples of each rejection reason
    print("\n" + "=" * 70)
    print("Examples by Filter Reason")
    print("=" * 70)

    from collections import defaultdict
    examples_by_reason = defaultdict(list)
    for result in batch_results:
        examples_by_reason[result.reason].append(result.pdb_id)

    for reason in FilterReason:
        pdbs = examples_by_reason.get(reason, [])
        if pdbs:
            print(f"\n{reason.value} ({len(pdbs)} structures):")
            for pdb in pdbs[:5]:  # Show first 5 examples
                print(f"  - {pdb}")
            if len(pdbs) > 5:
                print(f"  ... and {len(pdbs) - 5} more")

if __name__ == "__main__":
    main()
