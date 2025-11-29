# Module Dependencies

## Dependency Graph

```
                    ┌─────────────┐
                    │   cli.py    │
                    └──────┬──────┘
                           │
              ┌────────────┼────────────┐
              │            │            │
              ▼            ▼            ▼
        ┌──────────┐ ┌──────────┐ ┌───────────────┐
        │pipeline.py│ │yaml_conv.│ │ filtering.py  │
        └────┬─────┘ └──────────┘ └───────────────┘
             │
    ┌────────┼────────┬───────────┐
    │        │        │           │
    ▼        ▼        ▼           ▼
┌────────┐┌────────┐┌──────────┐┌──────────┐
│sabdab.py││mapping.││structure.││pymol_env.│
└────────┘└────────┘└──────────┘└──────────┘
```

## Detailed Dependencies

### sabdab.py
- **Purpose**: SAbDab data download and filtering
- **External deps**: pandas, requests, tqdm
- **Internal deps**: None (base module)
- **Used by**: pipeline.py

### mapping.py
- **Purpose**: UniProt/Ensembl API interactions, ortholog mapping
- **External deps**: requests, biopython (for sequence alignment)
- **Internal deps**: None (base module)
- **Used by**: pipeline.py

### structure.py
- **Purpose**: PDB/mmCIF download, parsing, alignment
- **External deps**: biopython, gemmi
- **Internal deps**: pymol_env.py (optional)
- **Used by**: pipeline.py

### pymol_env.py
- **Purpose**: PyMOL environment detection and setup
- **External deps**: None (subprocess for detection)
- **Internal deps**: None (base module)
- **Used by**: structure.py

### pipeline.py
- **Purpose**: Main orchestration, DataPoint management
- **External deps**: pandas, tqdm, gemmi
- **Internal deps**: sabdab.py, mapping.py, structure.py
- **Used by**: cli.py

### yaml_converter.py
- **Purpose**: Convert CIF files to Boltz YAML format
- **External deps**: gemmi, pyyaml
- **Internal deps**: None
- **Used by**: cli.py

### filtering.py
- **Purpose**: Filter data points by antibody-antigen contacts
- **External deps**: biopython
- **Internal deps**: None
- **Used by**: cli.py

### cli.py
- **Purpose**: Click-based command-line interface
- **External deps**: click
- **Internal deps**: pipeline.py, yaml_converter.py, filtering.py
- **Used by**: Entry point (antibody-abtigen command)

## External Dependencies Summary

| Package | Used By | Purpose |
|---------|---------|---------|
| biopython | mapping, structure, filtering | Sequence alignment, structure parsing |
| pandas | sabdab, pipeline | Data manipulation |
| requests | sabdab, mapping | HTTP API calls |
| tqdm | sabdab, pipeline | Progress bars |
| gemmi | structure, pipeline, yaml_converter | mmCIF parsing |
| pyyaml | yaml_converter | YAML generation |
| click | cli | CLI framework |

## Circular Dependencies

**Current status**: No circular dependencies detected.

The architecture follows a clean layered design:
1. Base modules (sabdab, mapping, structure, pymol_env) have no internal deps
2. Pipeline aggregates base modules
3. CLI consumes pipeline and utility modules
