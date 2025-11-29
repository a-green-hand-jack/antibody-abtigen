"""Pytest configuration and fixtures for antibody_abtigen tests."""
import pytest
from pathlib import Path
import shutil
import tempfile

# Test directories (isolated from main project)
TEST_DIR = Path(__file__).parent
TEST_DATA_DIR = TEST_DIR / "data"
TEST_OUTPUTS_DIR = TEST_DIR / "outputs"


@pytest.fixture(scope="session", autouse=True)
def setup_test_dirs():
    """Create test directories if they don't exist."""
    TEST_DATA_DIR.mkdir(exist_ok=True)
    TEST_OUTPUTS_DIR.mkdir(exist_ok=True)
    yield


@pytest.fixture
def test_data_dir():
    """Provide path to test data directory."""
    return TEST_DATA_DIR


@pytest.fixture
def test_outputs_dir():
    """Provide path to test outputs directory."""
    return TEST_OUTPUTS_DIR


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory for a single test."""
    temp_dir = TEST_OUTPUTS_DIR / "temp"
    temp_dir.mkdir(exist_ok=True)
    yield temp_dir
    # Cleanup after test
    if temp_dir.exists():
        shutil.rmtree(temp_dir)


@pytest.fixture
def temp_data_dir():
    """Create a temporary data directory for a single test."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_sabdab_row():
    """Provide a sample SAbDab row for testing."""
    import pandas as pd
    return pd.Series({
        'pdb': '5FUO',
        'Hchain': 'H',
        'Lchain': 'L',
        'antigen_chain': 'A',
        'antigen_species': 'homo sapiens',
        'antigen_type': 'protein',
        'antigen_name': 'Test Antigen',
        'resolution': 2.0,
    })


@pytest.fixture
def mock_uniprot_response():
    """Provide a mock UniProt API response."""
    return {
        'primaryAccession': 'P12345',
        'genes': [{'geneName': {'value': 'TEST_GENE'}}],
        'organism': {'taxonId': 9606, 'scientificName': 'Homo sapiens'},
        'sequence': {'value': 'MKTESTSEQUENCE'},
        'proteinDescription': {
            'recommendedName': {'fullName': {'value': 'Test Protein'}}
        },
        'uniProtKBCrossReferences': []
    }


@pytest.fixture
def mock_mouse_ortholog():
    """Provide a mock mouse ortholog response."""
    return {
        'accession': 'Q12345',
        'gene_name': 'Test',
        'taxon_id': 10090,
        'organism': 'Mus musculus',
        'sequence': 'MKTESTSEQUENCEMOUSE',
        'sequence_length': 19,
        'ensembl_ids': [],
        'protein_name': 'Test Protein Mouse',
        'sequence_identity': 85.0,
        'human_ortholog': 'P12345'
    }
