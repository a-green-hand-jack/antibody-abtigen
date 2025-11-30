"""
Test script for ESM-2 epitope encoder.

Tests:
1. Single structure encoding
2. Multi-chain epitope handling
3. HDF5 storage save/load
4. CSV logging
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from antibody_abtigen.epitope_pipeline import (
    GemmiStructureCleaner,
    GeometricEpitopeExtractor,
    ESM2EpitopeEncoder,
    HDF5EmbeddingStore,
    save_epitope_residues_csv,
    save_epitope_summary_csv,
    save_embedding_stats_csv,
    generate_epitope_report,
)


def main():
    print("=" * 70)
    print("ESM-2 Epitope Encoder Test")
    print("=" * 70)

    # Paths
    cleaned_cif_dir = Path("data/cleaned_cif")
    sabdab_summary = Path("data/meta/sabdab_summary_all.tsv")
    output_dir = Path("data/epitope_pipeline/embeddings")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get test CIF files (first 10, excluding any double-cleaned files)
    cif_files = [f for f in sorted(cleaned_cif_dir.glob("*_cleaned.cif"))
                 if "_CLEANED_" not in f.name][:10]
    print(f"\nTest files ({len(cif_files)}): {[f.stem for f in cif_files]}")

    # Initialize components
    print("\nInitializing components...")
    cleaner = GemmiStructureCleaner(
        sabdab_summary_path=sabdab_summary
    )
    extractor = GeometricEpitopeExtractor(distance_threshold=5.0)

    # Use cached model weights
    cache_dir = Path("/ibex/user/wuj0c/.cache")
    encoder = ESM2EpitopeEncoder(
        device="cuda",
        use_fp16=True,
        cache_dir=cache_dir
    )
    store = HDF5EmbeddingStore()

    # Process structures
    print("\nProcessing structures...")
    structures = {}
    epitopes = []
    encoder_outputs = []

    for cif_path in cif_files:
        # Extract PDB ID (will be normalized by cleaner)
        print(f"\n  Processing {cif_path.stem}...")

        try:
            # Clean structure (skip filtering since already cleaned)
            cleaned, filter_result = cleaner.clean_structure(cif_path, skip_filtering=True)
            if cleaned is None:
                print(f"    Skipped: {filter_result.reason.value}")
                continue

            # Use the normalized PDB ID from cleaner
            pdb_id = cleaned.pdb_id
            structures[pdb_id] = cleaned
            print(f"    PDB ID: {pdb_id}")
            print(f"    Chains: {[m.original_chain_id for m in cleaned.chain_mappings]}")

            # Extract epitope
            epitope = extractor.extract_epitope(cleaned)
            epitopes.append(epitope)
            print(f"    Epitope: {epitope.total_residue_count()} residues in {list(epitope.antigen_chains.keys())}")

            # Encode
            print(f"    Encoding with ESM-2...")
            output = encoder.encode_full(epitope, cleaned)
            encoder_outputs.append(output)

            print(f"    Full embedding shape: {output.full_embedding.shape}")
            print(f"    Epitope embedding shape: {output.epitope_embedding.shape}")
            print(f"    Chains encoded: {list(output.chain_embeddings.keys())}")

        except Exception as e:
            print(f"    Error: {e}")
            import traceback
            traceback.print_exc()
            continue

    if not encoder_outputs:
        print("\nNo structures encoded successfully!")
        return

    # Save embeddings to HDF5
    print("\n" + "=" * 70)
    print("Saving embeddings...")
    h5_path = output_dir / "test_embeddings.h5"
    store.save_encoder_outputs(encoder_outputs, h5_path)
    print(f"  Saved to: {h5_path}")

    # Get summary
    summary = store.get_summary(h5_path)
    print(f"  Summary: {summary}")

    # Save CSV logs
    print("\nSaving CSV logs...")

    csv_dir = output_dir / "csv"
    csv_dir.mkdir(exist_ok=True)

    save_epitope_residues_csv(
        epitopes, structures,
        csv_dir / "epitope_residues.csv"
    )
    print(f"  Saved: epitope_residues.csv")

    save_epitope_summary_csv(
        epitopes, structures,
        csv_dir / "epitope_summary.csv",
        encoder_outputs
    )
    print(f"  Saved: epitope_summary.csv")

    save_embedding_stats_csv(
        encoder_outputs,
        csv_dir / "embedding_stats.csv"
    )
    print(f"  Saved: embedding_stats.csv")

    # Generate report
    print("\n" + generate_epitope_report(epitopes, encoder_outputs))

    # Test loading
    print("\nTesting HDF5 load...")
    for output in encoder_outputs[:2]:
        loaded = store.load_encoder_output(h5_path, output.epitope_id)
        if loaded:
            print(f"  {output.epitope_id}: loaded successfully")
            print(f"    Full embedding norm: {loaded.full_embedding.sum():.4f}")
            print(f"    Epitope embedding norm: {loaded.epitope_embedding.sum():.4f}")
        else:
            print(f"  {output.epitope_id}: LOAD FAILED")

    print("\n" + "=" * 70)
    print("Test Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
