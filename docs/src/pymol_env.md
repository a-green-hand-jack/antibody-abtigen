# pymol_env.py - PyMOL Environment Module

## Overview

Detects and configures PyMOL for structure alignment across different platforms and installation methods.

## Location

`antibody_abtigen/pymol_env.py`

## Purpose

PyMOL provides better structure alignment quality than Biopython, but its installation varies by platform:
- macOS: PyMOL.app or conda
- Linux: conda/mamba only (pip package not available on most distros)
- Windows: pip or conda

This module handles detection and setup across all these scenarios.

## Public API

### Functions

#### `is_pymol_available() -> bool`

Check if PyMOL is available on the system.

#### `get_pymol_info() -> dict`

Get detailed PyMOL availability information.

**Returns:**
```python
{
    'available': True,
    'python_path': '/path/to/pymol/python',
    'method': 'macos_app'  # or 'conda', 'direct'
}
```

#### `setup_pymol(auto_create=False) -> Tuple[bool, str]`

Set up PyMOL environment.

**Parameters:**
- `auto_create`: If True, create conda environment if not found

**Returns:** Tuple of (success, message)

#### `get_pymol_python() -> Optional[str]`

Get path to PyMOL Python interpreter.

## Detection Order

1. **macOS PyMOL.app**:
   - `/Applications/PyMOL.app/Contents/bin/python`

2. **Direct Python import**:
   - `import pymol`

3. **Conda environments** (searched in order):
   - `pymol-env`
   - `pymol`
   - `pymol-open-source`
   - Any env with PyMOL installed

## Auto-Create Environment

On Linux/HPC systems, can auto-create conda environment:

```python
from antibody_abtigen.pymol_env import setup_pymol

success, message = setup_pymol(auto_create=True)
# Creates: mamba create -n pymol-env -c conda-forge pymol-open-source python=3.11 -y
```

## Usage Example

```python
from antibody_abtigen.pymol_env import (
    is_pymol_available,
    get_pymol_info,
    get_pymol_python
)

# Check availability
if is_pymol_available():
    info = get_pymol_info()
    print(f"PyMOL available via {info['method']}")
    print(f"Python: {info['python_path']}")
else:
    print("PyMOL not available, using Biopython fallback")

# Get Python path for subprocess calls
pymol_python = get_pymol_python()
if pymol_python:
    import subprocess
    subprocess.run([pymol_python, "-c", "import pymol; print('OK')"])
```

## HPC/Cluster Setup

For SLURM/HPC environments:

```bash
# One-time setup
mamba create -n pymol-env -c conda-forge pymol-open-source python=3.11 -y

# In job script
source activate pymol-env
# OR let the package auto-detect
uv run python run.py --limit 100
```

## Fallback Behavior

If PyMOL is unavailable:
- `structure.py` uses `align_structures_biopython()` instead
- Alignment quality may be lower but functional
- Use `--no-pymol` CLI flag to force Biopython
