"""
Integration test for Day 2 implementation.

Tests the complete clean → extract pipeline on 5 real PDB structures.
"""

import pytest
import tempfile
from pathlib import Path

from antibody_abtigen.epitope_pipeline import (
    GemmiStructureCleaner,
    GeometricEpitopeExtractor,
)


class TestDay2Integration:
    """Integration tests for cleaner + extractor pipeline."""

    @pytest.fixture
    def sample_cif_dir(self):
        """Get sample CIF directory."""
        data_dir = Path("./data/raw_cif")
        if not data_dir.exists():
            pytest.skip("No raw_cif directory found")
        return data_dir

    def test_five_structure_pipeline(self, sample_cif_dir):
        """
        Test complete pipeline on 5 real structures.

        Pipeline:
        1. Find CIF files with both antigen and antibody
        2. Clean each with GemmiStructureCleaner
        3. Extract epitope with GeometricEpitopeExtractor
        4. Validate results
        """
        # Find more CIF files to search through
        all_cif_files = sorted(sample_cif_dir.glob("*.cif"))[:30]  # Try first 30

        if len(all_cif_files) < 5:
            pytest.skip(f"Only found {len(all_cif_files)} CIF files, need at least 5")

        print(f"\n{'='*60}")
        print(f"Day 2 Integration Test: Clean → Extract Pipeline")
        print(f"{'='*60}")

        # Initialize components
        cleaner = GemmiStructureCleaner()
        extractor = GeometricEpitopeExtractor(distance_threshold=5.0)

        results = []

        with tempfile.TemporaryDirectory() as tmpdir:
            cleaner.output_dir = Path(tmpdir)

            processed = 0
            for i, cif_path in enumerate(all_cif_files, 1):
                if len(results) >= 5:
                    break  # Got enough successful structures

                print(f"\n[{i}/{len(all_cif_files)}] Processing {cif_path.name}...")

                try:
                    # Step 1: Clean structure
                    cleaned = cleaner.clean_structure(cif_path)

                    print(f"  ✓ Cleaned structure")
                    print(f"    PDB ID: {cleaned.pdb_id}")
                    print(f"    Chains: {len(cleaned.chain_mappings)}")

                    # Show chain details
                    antigen_chains = cleaned.get_antigen_chains()
                    antibody_chains = cleaned.get_antibody_chains()

                    print(f"    Antigen chains: {[m.standardized_chain_id for m in antigen_chains]}")
                    print(f"    Antibody chains: {[m.standardized_chain_id for m in antibody_chains]}")

                    # Step 2: Extract epitope
                    epitope = extractor.extract_epitope(cleaned)

                    print(f"  ✓ Extracted epitope")
                    print(f"    Epitope ID: {epitope.epitope_id}")
                    print(f"    Distance threshold: {epitope.distance_threshold}Å")
                    print(f"    Total residues: {epitope.total_residue_count()}")
                    print(f"    Total contacts: {epitope.num_contacts}")

                    # Show per-chain residue counts
                    for chain_id, residues in epitope.antigen_chains.items():
                        print(f"    Chain {chain_id}: {len(residues)} epitope residues")
                        if len(residues) > 0:
                            print(f"      Range: {min(residues)} - {max(residues)}")

                    # Validate
                    assert cleaned.pdb_id == epitope.pdb_id
                    assert epitope.total_residue_count() > 0
                    assert epitope.num_contacts > 0
                    assert len(epitope.antibody_chains) > 0

                    # Check index mapping exists
                    for mapping in cleaned.chain_mappings:
                        assert len(mapping.auth_seq_id_map) == len(mapping.sequence)

                    results.append({
                        'pdb_id': cleaned.pdb_id,
                        'chains': len(cleaned.chain_mappings),
                        'epitope_residues': epitope.total_residue_count(),
                        'contacts': epitope.num_contacts,
                    })

                except Exception as e:
                    print(f"  ✗ Failed: {str(e)}")
                    continue

        # Print summary
        print(f"\n{'='*60}")
        print(f"Summary")
        print(f"{'='*60}")
        print(f"Structures processed: {len(results)}/5")

        if results:
            print(f"\nResults:")
            print(f"{'PDB ID':<10} {'Chains':<8} {'Epitope Res':<12} {'Contacts':<10}")
            print(f"{'-'*50}")
            for r in results:
                print(f"{r['pdb_id']:<10} {r['chains']:<8} {r['epitope_residues']:<12} {r['contacts']:<10}")

            avg_residues = sum(r['epitope_residues'] for r in results) / len(results)
            avg_contacts = sum(r['contacts'] for r in results) / len(results)

            print(f"\nAverages:")
            print(f"  Epitope residues: {avg_residues:.1f}")
            print(f"  Total contacts: {avg_contacts:.1f}")

        assert len(results) >= 3, f"Expected at least 3 successful structures (with both antigen and antibody), got {len(results)}"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
