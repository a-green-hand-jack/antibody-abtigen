# Structure Cleaner Tests

This directory contains test scripts and outputs for the structure cleaner module.

## Test Scripts

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
