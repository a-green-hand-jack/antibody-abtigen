# structure.py - PDB Structure Operations Module

## Overview

Handles PDB/mmCIF structure downloading, parsing, chain extraction, and structural alignment using PyMOL or Biopython.

## Location

`antibody_abtigen/structure.py`

## Public API

### Constants

- `PYMOL_AVAILABLE`: Boolean indicating if PyMOL is available
- `PYMOL_PYTHON`: Path to PyMOL Python interpreter (if available)

### Classes

#### `ChainSelect(Select)`

Biopython Select subclass for filtering specific chains.

```python
io = PDBIO()
io.set_structure(structure)
io.save("output.pdb", ChainSelect(["A", "B"]))
```

### Functions

#### `download_structure(pdb_id: str, output_dir: str) -> Tuple[Optional[str], Optional[str]]`

Download PDB structure in both CIF and PDB formats.

**Parameters:**
- `pdb_id`: PDB ID (case-insensitive)
- `output_dir`: Directory to save files

**Returns:** Tuple of (cif_path, pdb_path), either may be None if download fails

**Caching:** Files are cached; existing files are not re-downloaded

#### `parse_structure(file_path: str, structure_id: str) -> Optional[Structure]`

Parse a PDB or mmCIF file into Biopython Structure object.

**Parameters:**
- `file_path`: Path to PDB or CIF file
- `structure_id`: ID to assign to structure

**Returns:** Biopython Structure or None

#### `save_structure(structure, output_path: str, chain_ids: List[str] = None)`

Save structure to PDB format, optionally filtering chains.

#### `save_structure_cif(structure, output_path: str, chain_ids: List[str] = None)`

Save structure to mmCIF format using Biopython.

#### `save_cif_with_metadata(input_cif, output_cif, chain_ids=None, rotran=None)`

Save CIF preserving original metadata, optionally with transformation.

**Parameters:**
- `input_cif`: Source mmCIF file
- `output_cif`: Output path
- `chain_ids`: Chains to keep (None = all)
- `rotran`: Tuple of (rotation_matrix, translation_vector) to apply

#### `align_structures_pymol(mobile_file, reference_file, output_file, mobile_chain, reference_chain) -> float`

Align structures using PyMOL's `super` command.

**Returns:** RMSD in Angstroms

**Note:** Requires PyMOL to be available; uses subprocess.

#### `align_structures_biopython(mobile_structure, reference_structure, mobile_chain, reference_chain) -> Tuple[float, Structure, Tuple]`

Align structures using Biopython Superimposer.

**Returns:** Tuple of (RMSD, aligned_structure, (rotation, translation))

#### `merge_structures(structure_chain_pairs, new_id) -> Structure`

Merge multiple structures/chains into a single structure.

**Parameters:**
- `structure_chain_pairs`: List of (structure, chain_ids) tuples
- `new_id`: ID for merged structure

## PyMOL Detection

The module automatically detects PyMOL from multiple sources:

1. **macOS PyMOL.app**: `/Applications/PyMOL.app/Contents/bin/python`
2. **Direct import**: `import pymol`
3. **Conda environments**: `pymol-env`, `pymol`, `pymol-open-source`

If PyMOL is unavailable, alignment falls back to Biopython.

## Usage Example

```python
from antibody_abtigen.structure import (
    download_structure,
    parse_structure,
    save_structure,
    align_structures_biopython,
    PYMOL_AVAILABLE
)

# Download structure
cif_path, pdb_path = download_structure("5FUO", "./cache")

# Parse structure
structure = parse_structure(cif_path, "5FUO")

# Save specific chains
save_structure(structure, "antibody.pdb", chain_ids=["H", "L"])

# Align two structures
rmsd, aligned, rotran = align_structures_biopython(
    mobile_structure=mouse_struct,
    reference_structure=human_struct,
    mobile_chain="A",
    reference_chain="A"
)
print(f"Alignment RMSD: {rmsd:.2f} A")
```
