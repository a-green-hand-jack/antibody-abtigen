"""
Unit tests for core data structures and configuration.
"""

import pytest
import numpy as np
from pathlib import Path
import tempfile

from antibody_abtigen.epitope_pipeline.core import (
    ChainMapping,
    CleanedStructure,
    EpitopeResidues,
    EpitopeEmbedding,
    EpitopeGroup,
    AlignedComplex,
    PipelineConfig,
    create_default_config,
    save_config_to_yaml,
    load_config_from_yaml,
    ConfigurationError
)


class TestChainMapping:
    """Test ChainMapping dataclass."""

    def test_create_chain_mapping(self):
        """Test basic creation."""
        mapping = ChainMapping(
            pdb_id="7K8T",
            original_chain_id="E",
            standardized_chain_id="A",
            chain_type="antigen",
            sequence="MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV",
            auth_seq_id_map={0: 100, 1: 101, 2: 102}
        )

        assert mapping.pdb_id == "7K8T"
        assert mapping.standardized_chain_id == "A"
        assert mapping.chain_type == "antigen"
        assert len(mapping.sequence) == 47  # Fixed: actual length
        assert mapping.auth_seq_id_map[0] == 100


class TestCleanedStructure:
    """Test CleanedStructure dataclass."""

    def test_get_antigen_chains(self):
        """Test filtering antigen chains."""
        mappings = [
            ChainMapping("7K8T", "E", "A", "antigen", "SEQ1"),
            ChainMapping("7K8T", "F", "B", "antigen", "SEQ2"),
            ChainMapping("7K8T", "H", "H", "antibody_heavy", "SEQ3"),
            ChainMapping("7K8T", "L", "L", "antibody_light", "SEQ4"),
        ]

        structure = CleanedStructure(
            pdb_id="7K8T",
            file_path=Path("/tmp/7K8T.cif"),
            chain_mappings=mappings,
            antigen_cluster_id="cluster_001"
        )

        antigen_chains = structure.get_antigen_chains()
        assert len(antigen_chains) == 2
        assert all(m.chain_type == "antigen" for m in antigen_chains)

    def test_get_antibody_chains(self):
        """Test filtering antibody chains."""
        mappings = [
            ChainMapping("7K8T", "E", "A", "antigen", "SEQ1"),
            ChainMapping("7K8T", "H", "H", "antibody_heavy", "SEQ2"),
            ChainMapping("7K8T", "L", "L", "antibody_light", "SEQ3"),
        ]

        structure = CleanedStructure(
            pdb_id="7K8T",
            file_path=Path("/tmp/7K8T.cif"),
            chain_mappings=mappings
        )

        ab_chains = structure.get_antibody_chains()
        assert len(ab_chains) == 2


class TestEpitopeResidues:
    """Test EpitopeResidues dataclass."""

    def test_single_chain_epitope(self):
        """Test epitope on single chain."""
        epitope = EpitopeResidues(
            epitope_id="7k8t_epi",
            pdb_id="7K8T",
            antigen_chains={"A": [100, 101, 102, 150, 151]},
            antibody_chains=["H", "L"],
            distance_threshold=5.0,
            num_contacts=25
        )

        assert epitope.total_residue_count() == 5
        residues = epitope.get_all_residues()
        assert len(residues) == 5
        assert residues[0] == ("A", 100)

    def test_multichain_epitope(self):
        """Test epitope spanning multiple chains."""
        epitope = EpitopeResidues(
            epitope_id="7k8t_epi",
            pdb_id="7K8T",
            antigen_chains={
                "A": [100, 101, 102],
                "B": [45, 46, 47]
            },
            antibody_chains=["H", "L"],
            num_contacts=30
        )

        assert epitope.total_residue_count() == 6
        residues = epitope.get_all_residues()
        assert ("A", 100) in residues
        assert ("B", 45) in residues


class TestEpitopeEmbedding:
    """Test EpitopeEmbedding dataclass."""

    def test_valid_embedding(self):
        """Test creating valid embedding."""
        # Create normalized vector
        vec = np.random.randn(2560).astype(np.float32)
        vec = vec / np.linalg.norm(vec)  # L2 normalize

        embedding = EpitopeEmbedding(
            epitope_id="7k8t_epi",
            pdb_id="7K8T",
            embedding=vec,
            full_sequence="MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDLLKRMKV",
            epitope_residue_indices=[0, 1, 2, 10, 11, 12]
        )

        assert embedding.embedding.shape == (2560,)
        assert np.isclose(np.linalg.norm(embedding.embedding), 1.0)

    def test_invalid_embedding_shape(self):
        """Test that wrong shape raises error."""
        vec = np.random.randn(512).astype(np.float32)  # Wrong size

        with pytest.raises(ValueError, match="must be shape"):
            EpitopeEmbedding(
                epitope_id="test",
                pdb_id="TEST",
                embedding=vec,
                full_sequence="SEQ",
                epitope_residue_indices=[0]
            )

    def test_invalid_embedding_normalization(self):
        """Test that unnormalized vector raises error."""
        vec = np.random.randn(2560).astype(np.float32)  # Not normalized

        with pytest.raises(ValueError, match="must be L2-normalized"):
            EpitopeEmbedding(
                epitope_id="test",
                pdb_id="TEST",
                embedding=vec,
                full_sequence="SEQ",
                epitope_residue_indices=[0]
            )


class TestEpitopeGroup:
    """Test EpitopeGroup dataclass."""

    def test_group_size(self):
        """Test group size calculation."""
        group = EpitopeGroup(
            group_id="group_001",
            reference_epitope_id="7k8t_epi",
            member_epitope_ids=["7k8t_epi", "8xyz_epi", "9abc_epi"],
            avg_similarity=0.92
        )

        assert group.size() == 3

    def test_get_mobile_ids(self):
        """Test getting IDs that need alignment."""
        group = EpitopeGroup(
            group_id="group_001",
            reference_epitope_id="7k8t_epi",
            member_epitope_ids=["7k8t_epi", "8xyz_epi", "9abc_epi"],
            avg_similarity=0.92
        )

        mobile = group.get_mobile_ids()
        assert len(mobile) == 2
        assert "7k8t_epi" not in mobile
        assert "8xyz_epi" in mobile


class TestPipelineConfig:
    """Test PipelineConfig validation and creation."""

    def test_valid_config(self):
        """Test creating valid config."""
        config = create_default_config(
            data_dir=Path("/tmp/test_data"),
            device="cuda"
        )

        errors = config.validate()
        assert len(errors) == 0
        assert config.device == "cuda"
        assert config.use_fp16 is True

    def test_invalid_threshold(self):
        """Test that invalid thresholds are caught."""
        config = create_default_config(Path("/tmp/test"))
        config.similarity_threshold = 1.5  # Invalid

        errors = config.validate()
        assert len(errors) > 0
        assert any("similarity_threshold" in e for e in errors)

    def test_invalid_device(self):
        """Test that invalid device is caught."""
        config = create_default_config(Path("/tmp/test"))
        config.device = "tpu"  # Invalid

        errors = config.validate()
        assert len(errors) > 0
        assert any("device" in e for e in errors)


class TestConfigSerialization:
    """Test YAML config saving and loading."""

    def test_round_trip(self):
        """Test saving and loading config."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create config
            original_config = create_default_config(
                data_dir=Path(tmpdir) / "data",
                device="cpu",
                similarity_threshold=0.90
            )

            # Save to YAML
            yaml_path = Path(tmpdir) / "config.yaml"
            save_config_to_yaml(original_config, yaml_path)

            assert yaml_path.exists()

            # Load back
            loaded_config = load_config_from_yaml(yaml_path)

            # Check key fields match
            assert loaded_config.device == "cpu"
            assert loaded_config.similarity_threshold == 0.90
            assert loaded_config.embedding_batch_size == 8

    def test_load_nonexistent_file(self):
        """Test loading non-existent file raises error."""
        with pytest.raises(ConfigurationError, match="not found"):
            load_config_from_yaml(Path("/nonexistent/config.yaml"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
