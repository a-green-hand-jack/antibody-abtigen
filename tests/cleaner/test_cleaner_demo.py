"""
Demo script to test structure cleaner and save output.

This script demonstrates the cleaning pipeline and saves results
to a persistent directory for inspection.
"""

from pathlib import Path
from antibody_abtigen.epitope_pipeline import GemmiStructureCleaner, GeometricEpitopeExtractor

def main():
    # Setup paths (relative to project root)
    raw_cif_dir = Path("./data/raw_cif")
    output_dir = Path("./tests/cleaner/output_cleaned")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find first few CIF files with antibody-antigen complexes
    cif_files = sorted(raw_cif_dir.glob("*.cif"))[:10]

    # Initialize cleaner and extractor
    sabdab_summary = Path("./data/meta/sabdab_summary_all.tsv")
    cleaner = GemmiStructureCleaner(
        output_dir=output_dir,
        sabdab_summary_path=sabdab_summary if sabdab_summary.exists() else None
    )
    extractor = GeometricEpitopeExtractor(distance_threshold=5.0)

    print("=" * 70)
    print("Structure Cleaner Demo")
    print("=" * 70)

    successful = []

    for i, cif_path in enumerate(cif_files, 1):
        print(f"\n[{i}/{len(cif_files)}] Processing {cif_path.name}...")

        try:
            # Clean structure
            cleaned = cleaner.clean_structure(cif_path)

            print(f"  ✓ Cleaned: {cleaned.file_path}")
            print(f"    PDB ID: {cleaned.pdb_id}")
            print(f"    Total chains: {len(cleaned.chain_mappings)}")

            # Show chain details
            antigen_chains = cleaned.get_antigen_chains()
            antibody_chains = cleaned.get_antibody_chains()

            print(f"    Antigen chains: {[m.standardized_chain_id for m in antigen_chains]}")
            print(f"    Antibody chains: {[m.standardized_chain_id for m in antibody_chains]}")

            # Show sequence lengths
            for mapping in cleaned.chain_mappings:
                print(f"      {mapping.standardized_chain_id} ({mapping.chain_type}): "
                      f"{len(mapping.sequence)} residues")

            # Try to extract epitope
            if antigen_chains and antibody_chains:
                epitope = extractor.extract_epitope(cleaned)
                print(f"    Epitope residues: {epitope.total_residue_count()}")
                successful.append(cleaned.pdb_id)
            else:
                print(f"    ⚠ Skipped epitope extraction (missing antigen or antibody)")

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            continue

    print("\n" + "=" * 70)
    print(f"Summary")
    print("=" * 70)
    print(f"Cleaned CIF files saved to: {output_dir}")
    print(f"Successful structures: {len(successful)}")

    if successful:
        print(f"\nYou can inspect these cleaned structures:")
        for pdb_id in successful[:5]:
            print(f"  - {output_dir}/{pdb_id}_cleaned.cif")

if __name__ == "__main__":
    main()
