# Structure Cleaner Tests

This directory contains test scripts and outputs for the structure cleaner module.

## Test Scripts

### `test_filtering.py`
Tests antigen type filtering functionality on edge cases (protein, protein+peptide, peptide-only, no antigen).

**Usage**:
```bash
uv run python tests/cleaner/test_filtering.py
```

**Features**:
- Validates filtering rules on specific test cases
- Tests batch filtering on first 50 structures
- Generates CSV filter log and summary statistics

### `batch_clean_with_filter.py`
Production script for batch cleaning with filtering on all CIF files.

**Usage**:
```bash
uv run python tests/cleaner/batch_clean_with_filter.py
```

**Features**:
- Processes all CIF files in `data/raw_cif/`
- Filters out non-protein antigens
- Saves only accepted structures to `data/epitope_pipeline/cleaned/`
- Generates comprehensive filtering statistics
- Extracts epitopes from accepted structures

### `test_cleaner_demo.py`
Main demo script that tests the complete clean → extract pipeline on 10 structures.

**Usage**:
```bash
uv run python tests/cleaner/test_cleaner_demo.py
```

**Features**:
- Uses SAbDab metadata for accurate chain classification
- Tests cleaning on first 10 CIF files
- Validates epitope extraction
- Saves cleaned structures to `output_cleaned/`

### `check_cif_chains.py`
Diagnostic script to compare chain IDs between original and cleaned CIF files.

**Usage**:
```bash
uv run python tests/cleaner/check_cif_chains.py
```

### `debug_cif_structure.py`
Low-level CIF structure debugging tool using Gemmi.

**Usage**:
```bash
uv run python tests/cleaner/debug_cif_structure.py
```

### `test_load_cleaned_cif.py`
Validates that cleaned CIF files can be properly loaded and have valid coordinates.

**Usage**:
```bash
uv run python tests/cleaner/test_load_cleaned_cif.py
```

### `convert_cleaned_to_pdb.py`
Converts cleaned CIF files to PDB format for easier visualization.

**Usage**:
```bash
uv run python tests/cleaner/convert_cleaned_to_pdb.py
```

## Output Directories

### `output_cleaned/`
Contains cleaned CIF files with:
- Original PDB chain IDs preserved (no renaming)
- Water molecules removed
- HETATM records removed
- Only antigen + antibody chains retained

### `output_pdb/`
Contains PDB format versions of cleaned structures for tools that prefer PDB over mmCIF.

## Filtering Logic

### Antigen Type Filtering

The cleaner now filters structures based on antigen type to ensure only suitable structures are processed for ESM-2 embedding:

**Accepted**:
- Protein antigens (pure protein)
- Protein + peptide combinations (protein component used for embedding)

**Rejected**:
1. No antigen (antigen_chain=NA) - No epitope to extract
2. Peptide-only antigens - Too short for reliable ESM-2 embeddings
3. Hapten-only antigens - Small molecules, not proteins
4. Insufficient residues - Less than 10 total antigen residues

**Filter Results**:
All filtering decisions are logged to CSV with:
- `pdb_id`: Structure identifier
- `accepted`: Yes/No
- `filter_reason`: Reason for acceptance/rejection
- `antigen_chains`: Chain IDs from SAbDab
- `antigen_type`: Type from SAbDab (protein/peptide/hapten)
- `num_antigen_chains`: Number of antigen chains found
- `total_antigen_residues`: Total residues in all antigen chains
- `num_antibody_chains`: Number of antibody chains
- `details`: Human-readable explanation

### Output Files

When running batch cleaning:
- **Cleaned CIF files**: `data/epitope_pipeline/cleaned/*.cif` (only accepted structures)
- **Filter log**: `data/epitope_pipeline/filtering_log.csv` (all structures with decisions)
- **Epitope stats**: `data/epitope_pipeline/epitope_stats.csv` (epitope info for accepted structures)

## Key Design Decisions

### No Chain ID Renaming
We **do NOT rename chain IDs** in cleaned structures. Original PDB chain IDs are preserved because:
1. Maintains structure integrity for PyMOL visualization
2. Easier to trace back to original PDB entries
3. Avoids potential issues with mmCIF label_asym_id vs auth_asym_id

### SAbDab Metadata Integration
Chain classification uses SAbDab summary TSV (`data/meta/sabdab_summary_all.tsv`) to accurately identify:
- `Hchain`: Antibody heavy chain
- `Lchain`: Antibody light chain
- `antigen_chain`: Antigen chains

Falls back to length-based heuristics only when SAbDab metadata is unavailable.

## Example Output

For structure `1A14`:
- **Original chains**: N (388 res), H (120 res), L (104 res)
- **SAbDab metadata**: Hchain=H, Lchain=L, antigen_chain=N
- **Classification**:
  - N → antigen (neuraminidase)
  - H → antibody_heavy
  - L → antibody_light
- **Cleaned file**: `output_cleaned/1A14_cleaned.cif`
- **Epitope residues**: 21 (extracted from chain N in contact with H+L)

## Verification

To verify cleaned structures are valid:

```bash
# Check all cleaned files load correctly
uv run python tests/cleaner/test_load_cleaned_cif.py

# Compare specific structure
uv run python tests/cleaner/check_cif_chains.py
```

Expected output: All structures should load with correct atom counts and valid coordinates.
