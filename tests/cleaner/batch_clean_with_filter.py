"""
Batch cleaning script with filtering.

This script processes all available CIF files and:
1. Filters out structures without protein antigens
2. Saves only accepted structures
3. Generates comprehensive filtering statistics
"""

from pathlib import Path
from antibody_abtigen.epitope_pipeline import (
    GemmiStructureCleaner,
    GeometricEpitopeExtractor,
    save_filter_log,
    generate_filter_summary,
)

def main():
    print("=" * 70)
    print("Batch Structure Cleaning with Filtering")
    print("=" * 70)

    raw_cif_dir = Path("data/raw_cif")
    output_dir = Path("data/epitope_pipeline/cleaned")
    output_dir.mkdir(parents=True, exist_ok=True)

    sabdab_summary = Path("data/meta/sabdab_summary_all.tsv")

    # Initialize cleaner and extractor
    cleaner = GemmiStructureCleaner(
        output_dir=output_dir,
        sabdab_summary_path=sabdab_summary
    )
    extractor = GeometricEpitopeExtractor(distance_threshold=5.0)

    # Get all CIF files
    cif_files = sorted(raw_cif_dir.glob("*.cif"))
    print(f"\nTotal CIF files found: {len(cif_files)}")

    filter_results = []
    cleaned_structures = []
    epitope_stats = []

    print("\nProcessing structures...")
    for i, cif_path in enumerate(cif_files, 1):
        if i % 100 == 0:
            print(f"  Processed {i}/{len(cif_files)}...")

        try:
            cleaned, filter_result = cleaner.clean_structure(cif_path)
            filter_results.append(filter_result)

            if cleaned:
                cleaned_structures.append(cleaned)

                # Try to extract epitope
                try:
                    epitope = extractor.extract_epitope(cleaned)
                    epitope_stats.append({
                        'pdb_id': cleaned.pdb_id,
                        'epitope_residues': epitope.total_residue_count(),
                        'num_contacts': epitope.num_contacts,
                        'antigen_chains': len(epitope.antigen_chains)
                    })
                except Exception as e:
                    print(f"  Warning: Failed to extract epitope for {cleaned.pdb_id}: {e}")

        except Exception as e:
            print(f"  Error processing {cif_path.name}: {e}")

    print(f"  Completed: {len(cif_files)}/{len(cif_files)}")

    # Save filter log
    log_path = Path("data/epitope_pipeline/filtering_log.csv")
    save_filter_log(filter_results, log_path)
    print(f"\n✓ Filter log saved to: {log_path}")

    # Save epitope statistics
    if epitope_stats:
        import csv
        epitope_log_path = Path("data/epitope_pipeline/epitope_stats.csv")
        epitope_log_path.parent.mkdir(parents=True, exist_ok=True)

        with open(epitope_log_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['pdb_id', 'epitope_residues', 'num_contacts', 'antigen_chains'])
            writer.writeheader()
            writer.writerows(epitope_stats)

        print(f"✓ Epitope statistics saved to: {epitope_log_path}")

    # Print summary
    print("\n" + generate_filter_summary(filter_results))

    # Print epitope statistics
    if epitope_stats:
        import statistics
        epitope_counts = [s['epitope_residues'] for s in epitope_stats]
        print("\n" + "=" * 70)
        print("Epitope Extraction Statistics (Accepted Structures)")
        print("=" * 70)
        print(f"Total structures with epitopes: {len(epitope_stats)}")
        print(f"Average epitope size: {statistics.mean(epitope_counts):.1f} residues")
        print(f"Median epitope size: {statistics.median(epitope_counts):.1f} residues")
        print(f"Min epitope size: {min(epitope_counts)} residues")
        print(f"Max epitope size: {max(epitope_counts)} residues")
        print("=" * 70)

    # Show cleaned structure locations
    print(f"\n✓ {len(cleaned_structures)} cleaned structures saved to:")
    print(f"  {output_dir}")
    if cleaned_structures:
        print(f"\nFirst 5 cleaned structures:")
        for cleaned in cleaned_structures[:5]:
            print(f"  - {cleaned.file_path.name}")

if __name__ == "__main__":
    main()
