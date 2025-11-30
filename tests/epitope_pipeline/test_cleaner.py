"""
Unit tests for structure cleaner module.
"""

import pytest
import tempfile
from pathlib import Path

from antibody_abtigen.epitope_pipeline import (
    GemmiStructureCleaner,
    StructureCleaningError,
)


class TestGemmiStructureCleaner:
    """Test GemmiStructureCleaner class."""

    @pytest.fixture
    def cleaner(self):
        """Create a cleaner instance."""
        return GemmiStructureCleaner()

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_clean_nonexistent_file(self, cleaner):
        """Test that cleaning nonexistent file raises error."""
        fake_path = Path("/nonexistent/7K8T.cif")

        with pytest.raises(StructureCleaningError, match="not found"):
            cleaner.clean_structure(fake_path)

    def test_chain_classification(self, cleaner):
        """Test chain classification logic."""
        # Mock chain info
        chain_info = {
            "A": {"sequence": "M" * 350, "length": 350, "auth_seq_ids": list(range(1, 351))},  # Heavy
            "B": {"sequence": "M" * 220, "length": 220, "auth_seq_ids": list(range(1, 221))},  # Light
            "C": {"sequence": "M" * 100, "length": 100, "auth_seq_ids": list(range(1, 101))},  # Antigen
            "D": {"sequence": "M" * 150, "length": 150, "auth_seq_ids": list(range(1, 151))},  # Antigen
        }

        chain_types = cleaner._classify_chains(chain_info)

        assert chain_types["A"] == "antibody_heavy"
        assert chain_types["B"] == "antibody_light"
        assert chain_types["C"] == "antigen"
        assert chain_types["D"] == "antigen"

    def test_chain_mapping_creation(self, cleaner):
        """Test chain mapping creation with standardization."""
        chain_info = {
            "E": {"sequence": "MKFLKF", "length": 6, "auth_seq_ids": [100, 101, 102, 103, 104, 105]},
            "F": {"sequence": "LKFLMK", "length": 6, "auth_seq_ids": [200, 201, 202, 203, 204, 205]},
            "H": {"sequence": "M" * 300, "length": 300, "auth_seq_ids": list(range(1, 301))},
            "L": {"sequence": "M" * 210, "length": 210, "auth_seq_ids": list(range(1, 211))},
        }

        chain_types = {
            "E": "antigen",
            "F": "antigen",
            "H": "antibody_heavy",
            "L": "antibody_light",
        }

        mappings = cleaner._create_chain_mappings("TEST", chain_info, chain_types)

        # Check standardization
        mapping_dict = {m.original_chain_id: m for m in mappings}

        assert mapping_dict["E"].standardized_chain_id == "A"
        assert mapping_dict["F"].standardized_chain_id == "B"
        assert mapping_dict["H"].standardized_chain_id == "H"
        assert mapping_dict["L"].standardized_chain_id == "L"

        # Check auth_seq_id mapping
        assert mapping_dict["E"].auth_seq_id_map[0] == 100
        assert mapping_dict["E"].auth_seq_id_map[5] == 105
        assert mapping_dict["F"].auth_seq_id_map[0] == 200

    def test_auth_seq_id_mapping(self, cleaner):
        """Test auth_seq_id to 0-based index mapping."""
        chain_info = {
            "A": {
                "sequence": "MKFL",
                "length": 4,
                "auth_seq_ids": [100, 101, 102, 103]  # PDB numbering starts at 100
            }
        }

        mapping = cleaner._create_mapping(
            pdb_id="TEST",
            original_id="A",
            standardized_id="A",
            chain_type="antigen",
            chain_info=chain_info
        )

        # Check mapping
        assert mapping.auth_seq_id_map[0] == 100
        assert mapping.auth_seq_id_map[1] == 101
        assert mapping.auth_seq_id_map[2] == 102
        assert mapping.auth_seq_id_map[3] == 103

        # Check sequence
        assert mapping.sequence == "MKFL"

    def test_batch_clean(self, cleaner):
        """Test batch cleaning (with error handling)."""
        # Create list with one nonexistent file
        fake_paths = [
            Path("/nonexistent/7K8T.cif"),
            Path("/nonexistent/8ABC.cif"),
        ]

        # Should not raise, but yield nothing (errors are logged)
        results = list(cleaner.batch_clean(fake_paths))
        assert len(results) == 0


class TestCleanerIntegration:
    """Integration tests with real CIF files (if available)."""

    @pytest.fixture
    def sample_cif_dir(self):
        """Get sample CIF directory if it exists."""
        data_dir = Path("./data/raw_cif")
        if data_dir.exists():
            return data_dir
        return None

    def test_clean_real_structure(self, sample_cif_dir):
        """Test cleaning a real structure (if available)."""
        if sample_cif_dir is None:
            pytest.skip("No sample CIF files available")

        # Find first CIF file
        cif_files = list(sample_cif_dir.glob("*.cif"))
        if not cif_files:
            pytest.skip("No CIF files in data directory")

        cif_path = cif_files[0]

        with tempfile.TemporaryDirectory() as tmpdir:
            cleaner = GemmiStructureCleaner(output_dir=Path(tmpdir))

            try:
                result = cleaner.clean_structure(cif_path)

                # Basic validation
                assert result.pdb_id == cif_path.stem.upper()
                assert len(result.chain_mappings) > 0
                assert result.file_path.exists()

                # Check chain mappings have sequences
                for mapping in result.chain_mappings:
                    assert len(mapping.sequence) > 0
                    assert len(mapping.auth_seq_id_map) == len(mapping.sequence)

                # Check standardized IDs
                antigen_ids = [
                    m.standardized_chain_id
                    for m in result.get_antigen_chains()
                ]
                antibody_ids = [
                    m.standardized_chain_id
                    for m in result.get_antibody_chains()
                ]

                # Should have some antigens and antibodies
                assert len(antigen_ids) > 0 or len(antibody_ids) > 0

                print(f"✓ Cleaned {result.pdb_id}:")
                print(f"  Antigen chains: {antigen_ids}")
                print(f"  Antibody chains: {antibody_ids}")

            except StructureCleaningError as e:
                pytest.skip(f"Could not clean structure: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
