# Epitope-Centric Pipeline Implementation Progress

**Date**: 2025-11-30
**Status**: Day 1 Complete ✅
**Branch**: `embedding`

## Overview

This document tracks the implementation of the epitope-centric antibody-antigen dataset pipeline as specified in `docs/dev/不考虑同源关系通过表位定义相似.md`.

The pipeline groups antibody-antigen structures by **epitope structural similarity** (using ESM-2 embeddings) rather than sequence homology or species relationships.

## Implementation Plan Summary

**Architecture**: Incremental extension approach
- New subpackage: `antibody_abtigen/epitope_pipeline/`
- New CLI command: `antibody-abtigen epitope-cluster`
- No changes to existing `build` command

**Timeline**: 8-10 working days

**Key Technologies**:
- ESM-2 3B model for epitope embeddings
- FAISS GPU for similarity search
- NetworkX for graph clustering
- PyMOL for epitope-based alignment
- Gemmi for CIF processing

## Day 1: Core Infrastructure (COMPLETED ✅)

### Implemented Modules

#### 1. `antibody_abtigen/epitope_pipeline/core/dataclasses.py` (204 lines)

**Data structures**:
- `ChainMapping`: Maps PDB chain IDs to standardized IDs with index mapping
- `CleanedStructure`: Result of structure cleaning (Phase 1)
- `EpitopeResidues`: Epitope definition with multi-chain support (Phase 2)
- `EpitopeEmbedding`: ESM-2 embedding with validation (Phase 3)
- `EpitopeGroup`: Similarity-based epitope grouping (Phase 4)
- `AlignedComplex`: PyMOL alignment result (Phase 5)
- `PipelineConfig`: Configuration with validation

**Key features**:
- Immutable dataclasses (frozen=True)
- Index mapping: 0-based (ESM-2) ↔ PDB auth_seq_id ↔ PyMOL resi
- L2 normalization validation for embeddings
- Configuration validation with detailed error messages

#### 2. `antibody_abtigen/epitope_pipeline/core/interfaces.py` (193 lines)

**Abstract base classes**:
- `StructureCleaner`: CIF cleaning and standardization
- `EpitopeExtractor`: Distance-based epitope extraction
- `EpitopeEncoder`: ESM-2 embedding generation
- `EmbeddingStore`: HDF5 storage interface
- `EpitopeGrouper`: FAISS + NetworkX clustering
- `StructureAligner`: PyMOL/Biopython alignment
- `PipelineOrchestrator`: High-level pipeline control

**Design principles**:
- Dependency injection for testability
- Easy swapping of implementations
- Clear separation of concerns

#### 3. `antibody_abtigen/epitope_pipeline/core/config.py` (227 lines)

**Configuration management**:
- YAML file loading/saving
- Nested structure flattening
- Default configuration generation
- Comprehensive validation

**Example config structure**:
```yaml
data:
  raw_cif_dir: data/raw_cif
  cleaned_cif_dir: data/epitope_pipeline/cleaned
  embeddings_dir: data/epitope_pipeline/embeddings

embedding:
  model_name: esm2_t36_3B_UR50D
  batch_size: 8
  device: cuda
  use_fp16: true

grouping:
  similarity_threshold: 0.85
  top_k_neighbors: 100

alignment:
  method: pymol
  use_super: true
```

#### 4. `antibody_abtigen/epitope_pipeline/core/exceptions.py` (34 lines)

**Custom exceptions**:
- `EpitopePipelineError` (base)
- `StructureCleaningError`
- `EpitopeExtractionError`
- `EmbeddingError`
- `GroupingError`
- `AlignmentError`
- `ConfigurationError`
- `IndexMappingError`

#### 5. `tests/epitope_pipeline/test_core.py` (257 lines)

**Test coverage**: 15 tests, all passing ✅
- Chain mapping creation
- Structure filtering (antigen/antibody chains)
- Epitope residue handling (single/multi-chain)
- Embedding validation (shape, normalization)
- Group operations
- Config validation (thresholds, device)
- YAML serialization round-trip

### Test Results

```bash
$ uv run pytest tests/epitope_pipeline/test_core.py -v

============================= test session starts ==============================
tests/epitope_pipeline/test_core.py::TestChainMapping::test_create_chain_mapping PASSED
tests/epitope_pipeline/test_core.py::TestCleanedStructure::test_get_antigen_chains PASSED
tests/epitope_pipeline/test_core.py::TestCleanedStructure::test_get_antibody_chains PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeResidues::test_single_chain_epitope PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeResidues::test_multichain_epitope PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeEmbedding::test_valid_embedding PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeEmbedding::test_invalid_embedding_shape PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeEmbedding::test_invalid_embedding_normalization PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeGroup::test_group_size PASSED
tests/epitope_pipeline/test_core.py::TestEpitopeGroup::test_get_mobile_ids PASSED
tests/epitope_pipeline/test_core.py::TestPipelineConfig::test_valid_config PASSED
tests/epitope_pipeline/test_core.py::TestPipelineConfig::test_invalid_threshold PASSED
tests/epitope_pipeline/test_core.py::TestPipelineConfig::test_invalid_device PASSED
tests/epitope_pipeline/test_core.py::TestConfigSerialization::test_round_trip PASSED
tests/epitope_pipeline/test_core.py::TestConfigSerialization::test_load_nonexistent_file PASSED

============================== 15 passed in 0.63s ==============================
```

### Public API

```python
from antibody_abtigen.epitope_pipeline import (
    # Data structures
    CleanedStructure,
    EpitopeResidues,
    EpitopeEmbedding,
    EpitopeGroup,
    AlignedComplex,
    PipelineConfig,

    # Configuration
    create_default_config,
    load_config_from_yaml,
    save_config_to_yaml,

    # Exceptions
    EpitopePipelineError,
    # ... etc
)

# Example usage
config = create_default_config(
    data_dir="./data",
    device="cuda",
    similarity_threshold=0.85
)
```

## Directory Structure

```
antibody_abtigen/
├── epitope_pipeline/
│   ├── __init__.py                   # Public API exports
│   ├── core/
│   │   ├── __init__.py
│   │   ├── dataclasses.py            # ✅ Implemented
│   │   ├── interfaces.py             # ✅ Implemented
│   │   ├── config.py                 # ✅ Implemented
│   │   └── exceptions.py             # ✅ Implemented
│   ├── cleaner.py                    # 🔜 Day 2
│   ├── extractor.py                  # 🔜 Day 2
│   ├── encoder.py                    # 🔜 Day 3-4
│   ├── storage.py                    # 🔜 Day 3-4
│   ├── grouper.py                    # 🔜 Day 5
│   ├── aligner.py                    # 🔜 Day 6-7
│   └── orchestrator.py               # 🔜 Day 8
└── cli.py                            # 🔜 Day 8 (extend)

tests/
└── epitope_pipeline/
    ├── __init__.py
    └── test_core.py                  # ✅ 15/15 tests passing
```

## Next Steps (Day 2)

### Module I: Structure Cleaner (`cleaner.py`)

**Goal**: Clean CIF files and standardize chain IDs

**Tasks**:
- Implement `GemmiStructureCleaner` class
- Reuse `structure.py::clean_structure()` logic
- Add chain ID renaming (antigen → A,B,...; antibody → H,L)
- Build auth_seq_id → 0-based index mapping
- Unit tests with small synthetic CIF files

**Success criteria**: Extract clean structures with proper index mappings

### Module II: Epitope Extractor (`extractor.py`)

**Goal**: Define epitope residues via geometric contacts

**Tasks**:
- Implement `GeometricEpitopeExtractor` class
- Reuse `epitope.py::get_epitope_residues()`
- Add multi-chain epitope merging
- Handle trimer interfaces
- Unit tests with known epitopes

**Success criteria**: Extract epitopes from CIF and save to JSON

### Integration Test

Run complete clean → extract pipeline on 5 real PDB structures.

## Key Design Decisions

### 1. Full-Context Embedding Strategy

**Decision**: Process entire antigen sequence through ESM-2, then extract epitope residue embeddings.

**Rationale**: Antibodies recognize epitopes in their folding context. Truncating to epitope-only loses structural information.

**Implementation**:
```python
# Full sequence through ESM-2
results = model(full_sequence_tokens)
# Extract epitope positions
epitope_vectors = results[epitope_indices]
# Mean pooling
embedding = torch.mean(epitope_vectors, dim=0)
```

### 2. Index Mapping Architecture

**Challenge**: Three different numbering systems:
- ESM-2: 0-based sequential indices
- PDB: auth_seq_id (can have gaps, insertions)
- PyMOL: uses auth_seq_id for selections

**Solution**: Store mapping in `ChainMapping.auth_seq_id_map`
```python
auth_seq_id_map: Dict[int, int]  # {0: 100, 1: 101, 2: 102, ...}
```

### 3. Multi-Chain Epitope Support

**Challenge**: Antibodies can bind at chain interfaces (e.g., trimers)

**Solution**: `EpitopeResidues.antigen_chains` as dict:
```python
antigen_chains = {
    "A": [100, 101, 102],  # Residues from chain A
    "B": [45, 46, 47]       # Residues from chain B
}
```

### 4. Validation-First Design

All data structures validate on creation:
- Embedding shape must be (2560,)
- Embedding must be L2-normalized
- Config thresholds must be in valid ranges
- Devices must be 'cuda' or 'cpu'

This catches errors early rather than failing deep in the pipeline.

## Environment

- **GPU**: NVIDIA A100-SXM4-80GB (CUDA 12.8)
- **PyMOL**: `/ibex/user/wuj0c/env/mamba/pymol-env`
- **ESM-2 Model**: Downloading to `/ibex/user/wuj0c/.cache/torch/hub/checkpoints/esm2_t36_3B_UR50D.pt`
- **Python**: 3.10.18 (uv environment)
- **All dependencies**: Already in `pyproject.toml`

## References

- **Design Document**: `docs/dev/不考虑同源关系通过表位定义相似.md`
- **Implementation Plan**: `/home/wuj0c/.claude/plans/goofy-questing-lake.md`
- **Existing Code to Reuse**:
  - `antibody_abtigen/epitope.py` (epitope extraction)
  - `antibody_abtigen/structure.py` (CIF processing, PyMOL)
  - `antibody_abtigen/pipeline.py` (orchestration patterns)

## Git History

- **2025-11-30**: Day 1 implementation - Core infrastructure complete
  - Added `epitope_pipeline/core/` module
  - 15 unit tests, all passing
  - Public API established
