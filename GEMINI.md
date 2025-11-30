# GEMINI.md - Project Context

## Project Overview

**Name**: `antibody-antigen` (Feature Branch: `feature-sabdab_summary`)
**Description**: Cross-species antibody-antigen structure dataset builder. This project aims to implement an epitope-centric SAbDab pipeline, utilizing ESM-2 for representation and FAISS for efficient grouping.
**Status**: Early Development / In Progress (as of Nov 30, 2025).

## Architecture & Goals

The project is designed to process antibody-antigen structural data from SAbDab. The planned pipeline includes:
1.  **Data Cleaning**: Parsing CIF/PDB files, renaming chains, and standardizing data.
2.  **Epitope Extraction**: Defining epitopes based on contact distance (< 4.5 Å).
3.  **Representation**: Using ESM-2 (3B parameter model) for sequence embeddings.
4.  **Grouping**: Clustering similar epitopes using FAISS and graph-based methods (NetworkX).
5.  **Alignment**: Structural alignment using PyMOL (or Biopython as fallback).

## Directory Structure

*   **`src/`**: Main source code directory.
    *   `antibody-antigen/`: (Currently empty) Intended location for the package modules.
*   **`scripts/`**: Helper scripts and reference implementations.
    *   `reference/summary.py`: A reference script for SAbDab data processing (adapted from DiffAb), handling chain mapping, residue masking, and PDB parsing.
    *   `download_sabdab.py`: Utility to download SAbDab summary files.
*   **`docs/`**: Project documentation.
    *   `dev/`: Developer notes and progress tracking (`implementation_progress.md`).
    *   `src/`: Planned module documentation (`pipeline.md`, `sabdab.md`, etc.).
*   **`data/`**: (Git-ignored) Local data storage for SAbDab summaries, raw CIF/PDB files, and processed outputs.

## Development Setup

### Prerequisites
*   Python >= 3.10
*   `uv` (suggested by `uv.lock`) or standard `pip`.
*   **Key Dependencies**: `biopython`, `pandas`, `torch`, `fair-esm`, `gemmi`, `faiss-gpu` (see `pyproject.toml` for full list).

### Installation
To set up the environment using `uv`:
```bash
uv sync
```
Or with pip:
```bash
pip install -e .[dev]
```

### Running Scripts
The project is currently utilizing scripts in `scripts/`.
*   **Download SAbDab Summary**:
    ```bash
    python scripts/download_sabdab.py
    ```
*   **Run Reference Summary Processing**:
    ```bash
    python scripts/reference/summary.py --summary_dir data/meta/sabdab_summary_all.tsv --output_dir data/processed/
    ```

## Conventions
*   **Code Style**: Adheres to standard Python practices (PEP 8).
*   **Testing**: `pytest` is configured in `pyproject.toml`.
*   **Formatting/Linting**: `ruff` is included in dev dependencies.
