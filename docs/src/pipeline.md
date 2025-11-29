# pipeline.py - Main Pipeline Module

## Overview

Orchestrates the complete cross-species dataset building workflow, coordinating SAbDab data, ortholog mapping, structure processing, and output generation.

## Location

`antibody_abtigen/pipeline.py`

## Public API

### Classes

#### `DataPoint`

Dataclass representing a single data point in the dataset.

**Fields:**
| Field | Type | Description |
|-------|------|-------------|
| `id` | str | Unique identifier (e.g., "DP_5FUO_A") |
| `pdb_id_human` | str | Human complex PDB ID |
| `pdb_id_mouse` | str | Mouse antigen PDB ID |
| `human_antigen_chain` | str | Primary human antigen chain |
| `mouse_antigen_chain` | str | Mouse antigen chain |
| `heavy_chain` | str | Antibody heavy chain |
| `light_chain` | str | Antibody light chain |
| `human_uniprot` | str | Human UniProt accession |
| `mouse_uniprot` | str | Mouse UniProt accession |
| `gene_name` | str | Gene name |
| `protein_name` | str | Protein name |
| `sequence_identity` | float | Human-mouse sequence identity (%) |
| `resolution_human` | float | Human structure resolution (Å) |
| `resolution_mouse` | float | Mouse structure resolution (Å) |
| `alignment_rmsd` | float | Structure alignment RMSD (Å) |
| `status` | str | "pending", "success", "failed", "dry_run" |
| `error_message` | str | Error description if failed |
| `created_at` | str | ISO timestamp |
| `human_antigen_chains` | List[str] | All antigen chains |

#### `CrossSpeciesDatasetPipeline`

Main pipeline class for building the dataset.

**Constructor:**
```python
pipeline = CrossSpeciesDatasetPipeline(
    data_dir="./data",           # Cache directory
    output_dir="./output",       # Output directory
    resolution_threshold=2.5,    # Max resolution (Å)
    sequence_identity_threshold=50.0,  # Min identity (%)
    use_pymol=True               # Use PyMOL for alignment
)
```

**Methods:**

##### `run(limit=None, dry_run=False) -> pd.DataFrame`

Run the complete pipeline.

**Parameters:**
- `limit`: Max entries to process (None = all)
- `dry_run`: If True, analyze only without downloading structures

**Returns:** DataFrame with all DataPoints

**Pipeline steps:**
1. Download and parse SAbDab data
2. Filter for human antigen complexes
3. For each complex:
   - Map antigen to UniProt
   - Find mouse ortholog
   - Check sequence identity threshold
   - Find PDB structure for mouse protein
   - Download and align structures (unless dry_run)
   - Save output files
4. Generate summary CSV and log

##### `log(message: str, level: str = "INFO")`

Log a message to processing log.

## Output Structure

For each successful data point, creates:

```
output/DP_XXXX_Y/
├── DP_XXXX_Y_antibody.pdb        # Antibody (H+L chains)
├── DP_XXXX_Y_antibody.cif
├── DP_XXXX_Y_human_ag.pdb        # Human antigen (primary chain)
├── DP_XXXX_Y_human_ag.cif
├── DP_XXXX_Y_human_ag_full.pdb   # Human antigen (all chains)
├── DP_XXXX_Y_human_ag_full.cif
├── DP_XXXX_Y_mouse_ag.pdb        # Aligned mouse antigen
├── DP_XXXX_Y_mouse_ag.cif
├── DP_XXXX_Y_mouse_ag_full.pdb   # Full aligned mouse structure
├── DP_XXXX_Y_mouse_ag_full.cif
├── DP_XXXX_Y_human_complex.pdb   # Antibody + human antigen
├── DP_XXXX_Y_human_complex.cif
├── DP_XXXX_Y_mouse_complex.pdb   # Antibody + aligned mouse antigen
├── DP_XXXX_Y_mouse_complex.cif
└── metadata.json                  # Data point metadata
```

## Usage Example

```python
from antibody_abtigen.pipeline import CrossSpeciesDatasetPipeline

# Create pipeline
pipeline = CrossSpeciesDatasetPipeline(
    data_dir="./data",
    output_dir="./output",
    resolution_threshold=2.5,
    sequence_identity_threshold=50.0,
    use_pymol=False  # Use Biopython
)

# Dry run first
dry_results = pipeline.run(limit=10, dry_run=True)
print(f"Found {len(dry_results)} candidates")

# Full run
results = pipeline.run(limit=5, dry_run=False)
successful = results[results['status'] == 'success']
print(f"Successfully processed {len(successful)} data points")
```

## Failure Modes

Common failure reasons tracked in `error_message`:

- "Could not map antigen chain to UniProt"
- "Could not retrieve UniProt info"
- "No mouse ortholog found"
- "Sequence identity X% below threshold Y%"
- "No PDB structure found for mouse ortholog"
- "Structure processing error: ..."
