# Cross-Species Antibody-Antigen Structure Dataset Builder

[![CI](https://github.com/a-green-hand-jack/antibody-abtigen/actions/workflows/ci.yml/badge.svg)](https://github.com/a-green-hand-jack/antibody-abtigen/actions/workflows/ci.yml)

A pipeline for building cross-species antibody-antigen structure datasets from SAbDab and PDB data. This tool creates triplets of `<Mouse Antigen> <Human Antigen> <Antibody Template>` structures for computational biology and machine learning applications.

## Features

- Downloads and parses SAbDab (Structural Antibody Database) data
- Maps human antigens to mouse orthologs via UniProt/Ensembl APIs
- Downloads and processes PDB/mmCIF structure files
- Aligns mouse antigen structures to human antigen positions using PyMOL or Biopython
- Generates organized output with PDB/CIF files and metadata

## Installation

Requires Python 3.10+ and [uv](https://docs.astral.sh/uv/).

### Option 1: Install from GitHub (recommended for using as a package)

```bash
# In your project, add as dependency
uv add git+https://github.com/a-green-hand-jack/antibody-abtigen.git

# Or with PyMOL support (macOS/Windows only, for better alignment)
uv add "antibody-abtigen[pymol] @ git+https://github.com/a-green-hand-jack/antibody-abtigen.git"
```

### Option 2: Clone and develop locally

```bash
# Clone the repository
git clone https://github.com/a-green-hand-jack/antibody-abtigen.git
cd antibody-abtigen

# Install dependencies
uv sync

# Or install with PyMOL support (macOS/Windows only)
uv sync --extra pymol
```

> **Note on PyMOL**: The `pymol-open-source` package from PyPI is only available for macOS, Windows, and some Linux distributions. On unsupported platforms (e.g., RHEL-based Linux), see the [PyMOL on Linux](#pymol-on-linux-hpc-clusters) section below.

### PyMOL on Linux (HPC clusters)

On Linux systems where `pymol-open-source` from PyPI is not available, the package can automatically detect and use PyMOL from a conda/mamba environment:

```bash
# Create a conda environment with PyMOL (one-time setup)
mamba create -n pymol-env -c conda-forge pymol-open-source python=3.11 -y

# The package will automatically detect and use this environment
```

You can also programmatically check and set up PyMOL:

```python
from antibody_abtigen import setup_pymol, is_pymol_available, get_pymol_info

# Check PyMOL status
print(get_pymol_info())
# {'available': True, 'python_path': '/path/to/pymol-env/bin/python', 'method': 'conda'}

# Or setup with auto-create (will create pymol-env if not found)
success, message = setup_pymol(auto_create=True)
print(message)
```

The package searches for PyMOL in this order:
1. macOS PyMOL.app (`/Applications/PyMOL.app`)
2. Direct Python import (`import pymol`)
3. Conda environments named `pymol-env`, `pymol`, or `pymol-open-source`
4. Any conda environment with PyMOL installed

## Usage

### As a CLI tool

```bash
# If installed as package
antibody-abtigen --help
antibody-abtigen --limit 10 --dry-run

# Or run directly from repo
uv run python run.py --limit 10 --dry-run
```

### As a Python library

```python
from antibody_abtigen import CrossSpeciesDatasetPipeline

# Create pipeline
pipeline = CrossSpeciesDatasetPipeline(
    data_dir="./data",
    output_dir="./output",
    resolution_threshold=2.5,
    sequence_identity_threshold=50.0,
    use_pymol=False  # Use Biopython for alignment
)

# Run with limit
result_df = pipeline.run(limit=10, dry_run=False)

# Check results
print(f"Successful: {len(result_df[result_df['status'] == 'success'])}")
```

### Using individual modules

```python
from antibody_abtigen import (
    download_sabdab_summary,
    filter_human_antigen_complexes,
    get_uniprot_from_pdb_chain,
    find_mouse_ortholog,
    download_structure,
    align_structures_biopython,
)

# Download SAbDab data
summary_path = download_sabdab_summary("./data")

# Get UniProt ID from PDB
uniprot_id = get_uniprot_from_pdb_chain("5FUO", "A")

# Find mouse ortholog
mouse_info = find_mouse_ortholog(uniprot_id)
print(f"Mouse ortholog: {mouse_info['accession']}")
print(f"Sequence identity: {mouse_info['sequence_identity']:.1f}%")
```

### CLI Examples

```bash
# Run with default settings (all data)
antibody-abtigen build

# Demo mode with limited entries
antibody-abtigen build --limit 10

# Dry-run mode (analysis only, no structure downloads)
antibody-abtigen build --limit 50 --dry-run

# Filter to retain only antibody–antigen contacts (after build)
antibody-abtigen filter-interactions --input ./output --output ./output_filtered

# Custom thresholds
antibody-abtigen build --resolution 3.0 --identity 60 --limit 100

# Use Biopython instead of PyMOL for alignment
antibody-abtigen build --no-pymol
```

### Epitope Pipeline (ESM-2 Embeddings + Structure Validation)

The epitope pipeline groups antibody-antigen complexes by epitope similarity using a 5-stage process:
1. **Clean**: Filter and process raw CIF structures
2. **Embed**: Generate ESM-2 embeddings for epitope residues
3. **Group**: Cluster epitopes by embedding similarity
4. **Validate**: Filter groups by structural RMSD of pocket regions (NEW)
5. **Align**: Align structures within validated groups

**Full Pipeline (Recommended)**:
```bash
# Run all 5 stages with structure validation
antibody-abtigen epitope-pipeline \
    --input ./data/raw_cif \
    --output ./data/epitope_output \
    --rmsd-threshold 3.0 \
    --min-coverage 0.5 \
    --limit 100
```

**Structure Validation Options**:
| Option | Description | Default |
|--------|-------------|---------|
| `--skip-validate` | Skip structure validation (use embedding-only groups) | False |
| `--rmsd-threshold` | Maximum pocket RMSD (Å) for structure validation | 3.0 |
| `--min-coverage` | Minimum fraction of pocket atoms that must align | 0.5 |
| `--neighborhood-size` | Pocket expansion window size (residues) | 10 |

The validation stage expands discrete epitope contact residues to contiguous pocket regions (inspired by Boltz-2 Affinity Cropper), then computes pairwise RMSD between pocket Cα atoms to ensure structural similarity within groups.

**Step-by-Step**:
```bash
# Step 1: Clean and filter raw CIF files
antibody-abtigen clean \
    --input ./data/raw_cif \
    --output ./data/cleaned_cif \
    --sabdab-summary ./data/meta/sabdab_summary_all.tsv

# Step 2: Generate ESM-2 embeddings for epitopes
antibody-abtigen embed \
    --input ./data/cleaned_cif \
    --output ./data/embeddings \
    --device cuda \
    --limit 100  # Optional: process first N structures

# Step 3: Group epitopes by embedding similarity
antibody-abtigen group \
    --input ./data/embeddings/embeddings.h5 \
    --output ./data/grouping \
    --threshold 0.85 \
    --save-matrix  # Optional: save full similarity matrix
```

Output from `embed` command:
- `embeddings.h5`: HDF5 file with full-chain and epitope embeddings (2560-dim, L2 normalized)
- `epitope_residues.csv`: Per-residue epitope info (chain_id, auth_seq_id, residue_type)
- `epitope_summary.csv`: Per-structure summary
- `embedding_stats.csv`: Embedding statistics

Output from `group` command:
- `groups.json`: Group members with detailed metadata (epitope residues, similarity scores)
- `grouping_stats.csv`: Per-group statistics (member count, avg/min/max similarity)
- `similarity_sparse.h5`: Pairwise similarities above threshold
- `similarity_matrix.h5`: Full similarity matrix (if `--save-matrix`)
- `grouping_report.txt`: Human-readable summary

**Note on validation**: When running the full pipeline, the validation stage will:
- Expand epitope residues to contiguous pocket regions
- Compute pairwise RMSD between pocket Cα atoms
- Filter out members exceeding the RMSD threshold
- Output `validated_groups.json` and `validation_report.csv` in the grouping directory

```bash
# Step 4: Align structures within groups
antibody-abtigen align \
    --groups ./data/grouping/groups.json \
    --structures ./data/cleaned_cif \
    --output ./data/aligned
```

Output from `align` command:
- `group_XXXX/reference/`: Reference structure (antigen + antibody CIF files)
- `group_XXXX/aligned/`: Aligned structures (transformed to reference position)
- `group_XXXX/group_metadata.json`: Alignment RMSD and transformation info
- `alignment_summary.csv`: Per-group alignment statistics

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--output, -o` | Output directory | `./output` |
| `--data-dir, -d` | Data cache directory | `./data` |
| `--limit, -n` | Max entries to process | All |
| `--resolution, -r` | Max resolution (Å) | 2.5 |
| `--identity, -i` | Min sequence identity (%) | 50.0 |
| `--dry-run` | Analysis only, no downloads | False |
| `--no-pymol` | Use Biopython for alignment | False |

## Output Structure

```
output/
├── dataset_summary.csv      # Summary of all data points
├── processing_log.json      # Detailed processing log
└── DP_XXXX_Y/               # Individual data point folders
    ├── DP_XXXX_Y_antibody.pdb
    ├── DP_XXXX_Y_antibody.cif
    ├── DP_XXXX_Y_human_ag.pdb
    ├── DP_XXXX_Y_human_ag.cif
    ├── DP_XXXX_Y_mouse_ag.pdb  # Aligned to human position
    ├── DP_XXXX_Y_mouse_ag.cif
    └── metadata.json
```

## Pipeline Overview

1. **SAbDab Download**: Fetches the SAbDab summary file containing antibody-antigen complex information
2. **Filtering**: Selects human antigen complexes with protein antigens and resolution d threshold
3. **Ortholog Mapping**: Uses UniProt/Ensembl APIs to find mouse orthologs for human antigens
4. **Structure Retrieval**: Downloads PDB structures for both human complexes and mouse antigens
5. **Alignment**: Superimposes mouse antigen onto human antigen position using sequence-based structural alignment
6. **Output Generation**: Creates organized folder structure with separate PDB/CIF files for antibody, human antigen, and aligned mouse antigen

## Requirements

- Python 3.10+
- Biopython
- pandas
- requests
- tqdm
- click
- PyMOL (optional, for better alignment quality)

## License

MIT License
