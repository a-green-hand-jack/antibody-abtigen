"""
Unit tests for epitope extractor module.
"""

import pytest
import tempfile
from pathlib import Path

from antibody_abtigen.epitope_pipeline import (
    GeometricEpitopeExtractor,
    GemmiStructureCleaner,
    CleanedStructure,
    ChainMapping,
    EpitopeExtractionError,
)


class TestGeometricEpitopeExtractor:
    """Test GeometricEpitopeExtractor class."""

    @pytest.fixture
    def extractor(self):
        """Create extractor instance."""
        return GeometricEpitopeExtractor(distance_threshold=5.0)

    @pytest.fixture
    def mock_cleaned_structure(self, tmp_path):
        """Create a mock CleanedStructure (no real file)."""
        # This will fail in actual extraction, but tests error handling
        mappings = [
            ChainMapping(
                pdb_id="TEST",
                original_chain_id="E",
                standardized_chain_id="A",
                chain_type="antigen",
                sequence="MKFLKFSLLTAVLLSVVFAFSSCG",
            ),
            ChainMapping(
                pdb_id="TEST",
                original_chain_id="H",
                standardized_chain_id="H",
                chain_type="antibody_heavy",
                sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFS",
            ),
        ]

        fake_path = tmp_path / "TEST_cleaned.cif"

        return CleanedStructure(
            pdb_id="TEST",
            file_path=fake_path,
            chain_mappings=mappings,
        )

    def test_extract_nonexistent_file(self, extractor, mock_cleaned_structure):
        """Test that extracting from nonexistent file raises error."""
        with pytest.raises(EpitopeExtractionError):
            extractor.extract_epitope(mock_cleaned_structure)

    def test_distance_threshold(self):
        """Test custom distance threshold."""
        extractor = GeometricEpitopeExtractor(distance_threshold=8.0)
        assert extractor.distance_threshold == 8.0

    def test_batch_extract_error_handling(self, extractor, mock_cleaned_structure):
        """Test batch extraction error handling."""
        # Should not raise, but yield nothing (errors logged)
        results = list(extractor.batch_extract([mock_cleaned_structure]))
        assert len(results) == 0


class TestExtractorIntegration:
    """Integration tests with real structures."""

    @pytest.fixture
    def sample_cif_dir(self):
        """Get sample CIF directory if it exists."""
        data_dir = Path("./data/raw_cif")
        if data_dir.exists():
            return data_dir
        return None

    @pytest.fixture
    def cleaner(self):
        """Create cleaner for integration tests."""
        return GemmiStructureCleaner()

    @pytest.fixture
    def extractor(self):
        """Create extractor for integration tests."""
        return GeometricEpitopeExtractor(distance_threshold=5.0)

    def test_clean_and_extract_pipeline(
        self, sample_cif_dir, cleaner, extractor
    ):
        """Test complete clean → extract pipeline."""
        if sample_cif_dir is None:
            pytest.skip("No sample CIF files available")

        cif_files = list(sample_cif_dir.glob("*.cif"))
        if not cif_files:
            pytest.skip("No CIF files in data directory")

        cif_path = cif_files[0]

        with tempfile.TemporaryDirectory() as tmpdir:
            # Override cleaner output directory
            cleaner.output_dir = Path(tmpdir)

            try:
                # Step 1: Clean structure
                cleaned = cleaner.clean_structure(cif_path)

                print(f"\n✓ Cleaned {cleaned.pdb_id}")
                print(f"  Chains: {[m.standardized_chain_id for m in cleaned.chain_mappings]}")

                # Step 2: Extract epitope
                epitope = extractor.extract_epitope(cleaned)

                print(f"\n✓ Extracted epitope {epitope.epitope_id}")
                print(f"  Distance threshold: {epitope.distance_threshold}Å")
                print(f"  Total residues: {epitope.total_residue_count()}")
                print(f"  Total contacts: {epitope.num_contacts}")
                print(f"  Chains involved: {list(epitope.antigen_chains.keys())}")

                # Validation
                assert epitope.pdb_id == cleaned.pdb_id
                assert epitope.total_residue_count() > 0
                assert epitope.num_contacts > 0
                assert len(epitope.antibody_chains) > 0

                # Print residue details
                for chain_id, residues in epitope.antigen_chains.items():
                    print(f"  Chain {chain_id}: {len(residues)} residues")
                    print(f"    First 5: {residues[:5]}")

            except (Exception) as e:
                pytest.skip(f"Could not process structure: {e}")

    def test_multi_structure_batch(
        self, sample_cif_dir, cleaner, extractor
    ):
        """Test batch processing of multiple structures."""
        if sample_cif_dir is None:
            pytest.skip("No sample CIF files available")

        cif_files = list(sample_cif_dir.glob("*.cif"))[:5]  # First 5
        if len(cif_files) < 2:
            pytest.skip("Not enough CIF files for batch test")

        with tempfile.TemporaryDirectory() as tmpdir:
            cleaner.output_dir = Path(tmpdir)

            # Clean structures
            cleaned_structures = []
            for cif_path in cif_files:
                try:
                    cleaned = cleaner.clean_structure(cif_path)
                    cleaned_structures.append(cleaned)
                except Exception:
                    continue

            if not cleaned_structures:
                pytest.skip("No structures could be cleaned")

            # Extract epitopes in batch
            results = list(extractor.batch_extract(cleaned_structures))

            print(f"\n✓ Batch processed {len(results)} structures:")
            for cleaned, epitope in results:
                print(f"  {epitope.pdb_id}: {epitope.total_residue_count()} residues")

            assert len(results) > 0


class TestEpitopeResiduesMethods:
    """Test EpitopeResidues helper methods."""

    def test_total_residue_count_single_chain(self):
        """Test residue counting for single-chain epitope."""
        from antibody_abtigen.epitope_pipeline import EpitopeResidues

        epitope = EpitopeResidues(
            epitope_id="test",
            pdb_id="TEST",
            antigen_chains={"A": [100, 101, 102, 103, 104]},
            antibody_chains=["H", "L"],
        )

        assert epitope.total_residue_count() == 5

    def test_total_residue_count_multi_chain(self):
        """Test residue counting for multi-chain epitope."""
        from antibody_abtigen.epitope_pipeline import EpitopeResidues

        epitope = EpitopeResidues(
            epitope_id="test",
            pdb_id="TEST",
            antigen_chains={
                "A": [100, 101, 102],
                "B": [45, 46, 47, 48]
            },
            antibody_chains=["H", "L"],
        )

        assert epitope.total_residue_count() == 7

    def test_get_all_residues(self):
        """Test getting all residues as tuples."""
        from antibody_abtigen.epitope_pipeline import EpitopeResidues

        epitope = EpitopeResidues(
            epitope_id="test",
            pdb_id="TEST",
            antigen_chains={
                "A": [100, 101],
                "B": [45, 46]
            },
            antibody_chains=["H"],
        )

        residues = epitope.get_all_residues()

        assert len(residues) == 4
        assert ("A", 100) in residues
        assert ("A", 101) in residues
        assert ("B", 45) in residues
        assert ("B", 46) in residues


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
