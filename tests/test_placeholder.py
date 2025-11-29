"""Placeholder test to verify pytest setup works."""


def test_import_package():
    """Test that the package can be imported."""
    import antibody_abtigen
    assert antibody_abtigen is not None


def test_import_pipeline():
    """Test that pipeline module can be imported."""
    from antibody_abtigen.pipeline import CrossSpeciesDatasetPipeline
    assert CrossSpeciesDatasetPipeline is not None


def test_import_sabdab():
    """Test that sabdab module can be imported."""
    from antibody_abtigen.sabdab import (
        download_sabdab_summary,
        parse_sabdab_summary,
        filter_human_antigen_complexes,
    )
    assert download_sabdab_summary is not None
    assert parse_sabdab_summary is not None
    assert filter_human_antigen_complexes is not None


def test_import_mapping():
    """Test that mapping module can be imported."""
    from antibody_abtigen.mapping import (
        get_uniprot_from_pdb_chain,
        find_mouse_ortholog,
    )
    assert get_uniprot_from_pdb_chain is not None
    assert find_mouse_ortholog is not None


def test_import_structure():
    """Test that structure module can be imported."""
    from antibody_abtigen.structure import (
        download_structure,
        parse_structure,
        align_structures_biopython,
    )
    assert download_structure is not None
    assert parse_structure is not None
    assert align_structures_biopython is not None


def test_test_directories_exist(test_data_dir, test_outputs_dir):
    """Test that isolated test directories are available."""
    assert test_data_dir.exists()
    assert test_outputs_dir.exists()
