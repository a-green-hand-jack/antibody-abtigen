# Implementation Progress Tracker

**Goal**: Implement the Epitope-Centric SAbDab Pipeline (A100 Optimized) based on `docs/dev/不考虑同源关系通过表位定义相似.md`.

## Status Overview
- **Start Date**: 2025-11-30
- **Current Status**: In Progress

## Modules Checklist

### 1. Environment & Setup
- [ ] Update `pyproject.toml` with dependencies (`torch`, `fair-esm`, `h5py`, `networkx`, `faiss-gpu`, `biotite`, `typer`).
- [ ] Create directory structure (`data/`, `models/`).
- [ ] Create scripts for background downloading (SAbDab, ESM-2 weights).

### 2. Data Cleaning (`cleaner.py`)
- [ ] Implement CIF reading and cleaning with `gemmi`.
- [ ] Implement chain renaming (Antigen->A, Antibody->H/L).
- [ ] Implement antigen clustering and metadata generation.

### 3. Epitope Extraction (`extractor.py`)
- [ ] Implement contact-based epitope definition (< 4.5 Å).
- [ ] Generate `epitope_def.json`.

### 4. Representation (`encoder.py`)
- [ ] Implement ESM-2 (3B) loading.
- [ ] Implement Full Context Embedding + Mean Pooling.
- [ ] Support Batch processing on GPU.

### 5. Grouping (`grouper.py`)
- [ ] Implement FAISS indexing and retrieval.
- [ ] Implement NetworkX graph construction and connected components.

### 6. Alignment (`aligner.py` & `scripts/run_pymol_align.py`)
- [ ] Implement `scripts/run_pymol_align.py` (PyMOL-env standalone script).
- [ ] Implement `aligner.py` wrapper in main environment.

### 7. Integration (`main.py`)
- [ ] Update CLI to support the full pipeline.

## Log
- **2025-11-30**: Initialized progress tracker.
