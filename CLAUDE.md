# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics pipeline for building cross-species antibody-antigen structure datasets. It creates triplets of `<Mouse Antigen> <Human Antigen> <Antibody Template>` structures by:

1. Downloading antibody-antigen complex data from SAbDab
2. Mapping human antigens to mouse orthologs via UniProt/Ensembl APIs
3. Downloading and aligning PDB structures
4. Generating organized datasets with metadata

The project is both a CLI tool (`antibody-abtigen`) and a Python library.

## Development Commands

### Environment Setup

```bash
# Install dependencies
uv sync

# Install with PyMOL support (macOS/Windows)
uv sync --extra pymol

# Install dev dependencies
uv sync --all-extras
```

### Running the Pipeline

```bash
# Run via convenience script (defaults to 'build' subcommand)
uv run python run.py --limit 10 --dry-run

# Run via CLI directly
uv run antibody-abtigen build --limit 10 --dry-run

# Convert CIF to YAML format
uv run antibody-abtigen to-yaml --input ./output --output ./output_yamls

# Filter for antibody-antigen contacts
uv run antibody-abtigen filter-interactions --input ./output --output ./output_filtered

# Epitope pipeline: full pipeline (recommended)
uv run antibody-abtigen epitope-pipeline --input ./data/raw_cif --output ./data/epitope_output --limit 10

# Or step-by-step:
uv run antibody-abtigen clean --input ./data/raw_cif --output ./data/cleaned_cif
uv run antibody-abtigen embed --input ./data/cleaned_cif --output ./data/embeddings --limit 10
uv run antibody-abtigen group --input ./data/embeddings/embeddings.h5 --output ./data/grouping --threshold 0.85
uv run antibody-abtigen align --groups ./data/grouping/groups.json --structures ./data/cleaned_cif --output ./data/aligned
```

### Testing

```bash
# Run basic dry-run test (no structure downloads)
uv run python run.py --limit 10 --dry-run --no-pymol

# Run full test with structure downloads (5 samples)
uv run python run.py --limit 5 --no-pymol

# Test interaction filtering
uv run antibody-abtigen filter-interactions --input ./output --output ./output_filtered --dry-run
```

### Linting

```bash
# Run ruff linter
uv run ruff check antibody_abtigen/ run.py --ignore E501
```

### CI/CD

The project uses GitHub Actions with three jobs:
- **test**: Dry-run test with 10 samples (runs on all PRs)
- **test-with-structures**: Full test with 5 samples + interaction filter (main branch only)
- **lint**: Code quality checks with ruff

## Code Architecture

### Module Structure

The `antibody_abtigen/` package contains:

- **`sabdab.py`**: SAbDab data download and filtering
  - Downloads summary TSV from Oxford's SAbDab database
  - Filters for human antigen complexes with protein antigens

- **`mapping.py`**: UniProt/Ensembl API interactions
  - Maps PDB chains to UniProt IDs via SIFTS
  - Finds mouse orthologs via Ensembl REST API
  - Calculates sequence identity between orthologs
  - Implements rate limiting for API calls

- **`structure.py`**: PDB/mmCIF structure operations
  - Downloads structures from RCSB PDB
  - Parses and extracts chains using Biopython
  - Performs structure alignment (PyMOL or Biopython)
  - Saves structures in PDB and mmCIF formats

- **`pymol_env.py`**: PyMOL environment detection
  - Detects PyMOL from multiple sources (macOS app, conda, direct import)
  - Auto-creates `pymol-env` conda environment if needed
  - Search priority: macOS PyMOL.app → direct import → conda envs

- **`pipeline.py`**: Main orchestration logic
  - `CrossSpeciesDatasetPipeline` class coordinates all steps
  - `DataPoint` dataclass represents individual structure triplets
  - Generates summary CSV and processing logs

- **`yaml_converter.py`**: CIF to Boltz YAML conversion
  - Extracts sequences and disulfide bonds from mmCIF files
  - Generates YAML configs for Boltz-1 structure prediction

- **`epitope_pipeline/`**: Epitope-centric dataset pipeline (ESM-2 embeddings)
  - `cleaner.py`: Structure cleaning with SAbDab metadata-based chain classification
  - `extractor.py`: Distance-based epitope extraction (default: 5.0Å)
  - `encoder.py`: ESM-2 3B epitope encoder with multi-chain support
  - `storage.py`: HDF5 embedding storage with gzip compression
  - `epitope_log.py`: CSV logging utilities for epitope residues

- **`filtering.py`**: Post-processing interaction filter
  - Validates antibody-antigen contacts via distance calculations
  - Removes data points without sufficient interactions

- **`cli.py`**: Click-based command-line interface
  - Three subcommands: `build`, `to-yaml`, `filter-interactions`
  - Entry point: `antibody_abtigen.cli:main`

### Data Flow

1. `sabdab.py` → Download and filter SAbDab summary
2. `mapping.py` → For each human antigen chain:
   - Get UniProt ID from PDB chain
   - Find mouse ortholog via Ensembl
   - Find PDB structures for mouse protein
3. `structure.py` → Download human complex and mouse antigen structures
4. `structure.py` → Align mouse antigen to human antigen position
5. `pipeline.py` → Save triplet (antibody, human_ag, aligned_mouse_ag) with metadata
6. (Optional) `yaml_converter.py` → Convert to Boltz YAML format
7. (Optional) `filtering.py` → Filter for actual contacts

### Output Structure

```
output/
├── dataset_summary.csv       # All data points with metadata
├── processing_log.json        # Detailed processing log
└── DP_XXXX_Y/                 # One folder per data point
    ├── DP_XXXX_Y_antibody.pdb
    ├── DP_XXXX_Y_antibody.cif
    ├── DP_XXXX_Y_human_ag.pdb
    ├── DP_XXXX_Y_human_ag.cif
    ├── DP_XXXX_Y_mouse_ag.pdb  # Aligned to human position
    ├── DP_XXXX_Y_mouse_ag.cif
    └── metadata.json
```

## PyMOL Handling

PyMOL is optional for better alignment quality. The codebase has sophisticated detection logic:

1. **Detection order** (in `structure.py` and `pymol_env.py`):
   - macOS: `/Applications/PyMOL.app/Contents/bin/python`
   - Direct import: `import pymol`
   - Conda environments: `pymol-env`, `pymol`, `pymol-open-source`

2. **Automatic fallback**: If PyMOL unavailable, uses Biopython for alignment

3. **HPC/Linux setup**: On clusters where PyMOL.app doesn't exist, the code can auto-create a conda environment:
   ```bash
   mamba create -n pymol-env -c conda-forge pymol-open-source python=3.11 -y
   ```

4. When working with PyMOL code, be aware:
   - Alignment uses subprocess calls to external PyMOL Python
   - RMSD is returned from alignment operations
   - Temporary files are used for PyMOL scripts

## API Rate Limiting

All external API calls are rate-limited via `RateLimiter` class in `mapping.py`:
- UniProt: 3 requests/second
- Ensembl: 10 requests/second
- SIFTS: 5 requests/second

Results are cached using `@lru_cache` decorators to minimize API calls.

## Common Development Patterns

### Running a Quick Test

```bash
# Always use --dry-run first to validate logic without downloads
uv run python run.py --limit 10 --dry-run --no-pymol

# Then run a small batch with actual downloads
uv run python run.py --limit 5 --no-pymol
```

### Adding a New Feature

1. Identify which module(s) to modify based on the architecture above
2. Add unit tests if modifying core logic in `mapping.py` or `structure.py`
3. Test with `--limit 5` before running on full dataset
4. Update CLI in `cli.py` if adding user-facing options
5. Update the main README.md if changing user interface

### Debugging Pipeline Failures

1. Check `output/processing_log.json` for detailed error traces
2. Check `output/dataset_summary.csv` for status column values
3. Use `--limit 1` to isolate a single failing case
4. Add debug prints in the relevant module (sabdab/mapping/structure)

## HPC/SLURM Usage

For running on clusters, see `scripts/run_full_pipeline.sh`:

```bash
# Submit to SLURM
sbatch scripts/run_full_pipeline.sh

# Or run interactively
bash scripts/run_full_pipeline.sh
```

The script handles:
- Virtual environment activation
- PyMOL environment detection/creation
- Two-step pipeline: build → convert to YAML
- Summary statistics generation

## Important Conventions

- **Always use `uv run`** for Python script execution (per global CLAUDE.md instructions)
- Structure alignment returns RMSD in Angstroms
- PDB IDs are case-insensitive but stored uppercase
- Chain IDs are case-sensitive
- Sequence identity is stored as percentage (0-100), not fraction
- The `run.py` convenience script defaults to the `build` subcommand if no subcommand is provided
