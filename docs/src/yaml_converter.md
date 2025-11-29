# yaml_converter.py - CIF to YAML Conversion Module

## Overview

Converts mmCIF structure files to Boltz-1 YAML configuration format for structure prediction.

## Location

`antibody_abtigen/yaml_converter.py`

## Purpose

Boltz-1 is a protein structure prediction tool that requires YAML configuration files describing sequences and structural constraints. This module extracts this information from mmCIF files.

## Public API

### Functions

#### `batch_convert_to_yamls(raw_data_dir, output_dir, include_bonds=True, verbose=True) -> dict`

Convert all CIF files in DP_* folders to YAML format.

**Parameters:**
- `raw_data_dir`: Directory containing DP_* folders
- `output_dir`: Directory for YAML output
- `include_bonds`: Include disulfide bond information
- `verbose`: Print progress

**Returns:** Dictionary with conversion statistics

#### `cif_to_yaml(cif_path, output_path, include_bonds=True)`

Convert a single CIF file to YAML format.

**Parameters:**
- `cif_path`: Path to input mmCIF file
- `output_path`: Path for output YAML file
- `include_bonds`: Include disulfide bonds

## YAML Output Format

```yaml
# Example output for Boltz-1
version: 1
sequences:
  - protein:
      id: A
      sequence: MKTESTSEQUENCE...
  - protein:
      id: B
      sequence: ANOTHERTESTSEQ...
constraints:
  bonds:
    - source: [A, 22]
      target: [A, 96]
      type: disulfide
```

## Extracted Information

From mmCIF files:
- **Sequences**: Amino acid sequences for each chain
- **Chain IDs**: Label and auth chain identifiers
- **Disulfide bonds**: Cysteine-cysteine bonds from `_struct_conn` category

## Dependencies

- `gemmi`: For mmCIF parsing
- `pyyaml`: For YAML generation

## Usage Example

```python
from antibody_abtigen.yaml_converter import batch_convert_to_yamls, cif_to_yaml

# Convert single file
cif_to_yaml(
    "output/DP_5FUO_A/DP_5FUO_A_antibody.cif",
    "yamls/DP_5FUO_A_antibody.yaml"
)

# Batch convert
results = batch_convert_to_yamls(
    raw_data_dir="./output",
    output_dir="./yamls",
    include_bonds=True
)
print(f"Converted: {results['success']}, Failed: {results['failed']}")
```

## CLI Usage

```bash
antibody-abtigen to-yaml \
  --input ./output \
  --output ./yamls \
  # --no-bonds  # Optional: exclude disulfide bonds
```

## Output Structure

```
yamls/
├── DP_5FUO_A_antibody.yaml
├── DP_5FUO_A_human_ag.yaml
├── DP_5FUO_A_mouse_ag.yaml
├── DP_5FUO_A_human_complex.yaml
├── DP_5FUO_A_mouse_complex.yaml
└── conversion_summary.json
```
