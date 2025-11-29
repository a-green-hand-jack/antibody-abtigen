# mapping.py - UniProt/Ensembl Mapping Module

## Overview

Handles mapping between PDB chains, UniProt accessions, and cross-species ortholog discovery via external APIs.

## Location

`antibody_abtigen/mapping.py`

## Public API

### Classes

#### `RateLimiter`

Simple rate limiter to avoid API throttling.

```python
limiter = RateLimiter(requests_per_second=3.0)
limiter.wait()  # Blocks if needed to maintain rate limit
```

### Functions

#### `get_uniprot_from_pdb_chain(pdb_id: str, chain_id: str) -> Optional[str]`

Get UniProt accession from PDB ID and chain ID using SIFTS.

**Parameters:**
- `pdb_id`: PDB ID (e.g., "5FUO")
- `chain_id`: Chain ID (e.g., "A")

**Returns:** UniProt accession or None

**Caching:** Results cached via `@lru_cache(maxsize=1000)`

#### `get_uniprot_info(accession: str) -> Optional[Dict]`

Get detailed UniProt entry information.

**Returns dict with:**
- `accession`: UniProt ID
- `gene_name`: Primary gene name
- `taxon_id`: NCBI taxonomy ID
- `organism`: Scientific name
- `sequence`: Amino acid sequence
- `sequence_length`: Length of sequence
- `ensembl_ids`: Cross-referenced Ensembl IDs
- `protein_name`: Recommended protein name

#### `find_mouse_ortholog(human_uniprot_id: str) -> Optional[Dict]`

Find mouse ortholog for a human UniProt entry.

**Strategy:**
1. Get human protein info and gene name
2. Search UniProt for mouse proteins with same gene name
3. Calculate sequence identity between human and mouse
4. Return best match with identity score

**Returns dict with:**
- All fields from `get_uniprot_info()`
- `sequence_identity`: Percentage identity (0-100)
- `human_ortholog`: Original human UniProt ID

#### `find_pdb_structures_for_uniprot(uniprot_id: str) -> List[Dict]`

Find PDB structures containing a given UniProt protein.

**Returns list of dicts with:**
- `pdb_id`: PDB ID
- `entity_id`: Entity identifier
- `uniprot_id`: UniProt accession

#### `get_pdb_chain_for_entity(pdb_id: str, entity_id: str) -> List[str]`

Get chain IDs for a given PDB entity.

#### `calculate_sequence_identity(seq1: str, seq2: str) -> float`

Calculate sequence identity using Biopython global alignment.

**Returns:** Percentage identity (0-100)

## Constants

- `HUMAN_TAXON = 9606`: NCBI taxonomy ID for Homo sapiens
- `MOUSE_TAXON = 10090`: NCBI taxonomy ID for Mus musculus

## API Endpoints

| API | URL Pattern | Rate Limit |
|-----|-------------|------------|
| SIFTS | `ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}` | 5 req/s |
| UniProt | `rest.uniprot.org/uniprotkb/...` | 3 req/s |
| Ensembl | `rest.ensembl.org/homology/id/...` | 10 req/s |
| RCSB PDB | `search.rcsb.org/rcsbsearch/v2/query` | - |

## Usage Example

```python
from antibody_abtigen.mapping import (
    get_uniprot_from_pdb_chain,
    get_uniprot_info,
    find_mouse_ortholog,
    find_pdb_structures_for_uniprot
)

# Map PDB chain to UniProt
uniprot_id = get_uniprot_from_pdb_chain("5FUO", "A")
print(f"UniProt: {uniprot_id}")

# Get protein info
info = get_uniprot_info(uniprot_id)
print(f"Gene: {info['gene_name']}")

# Find mouse ortholog
mouse = find_mouse_ortholog(uniprot_id)
print(f"Mouse ortholog: {mouse['accession']}")
print(f"Sequence identity: {mouse['sequence_identity']:.1f}%")

# Find PDB structures for mouse protein
structures = find_pdb_structures_for_uniprot(mouse['accession'])
print(f"Found {len(structures)} PDB structures")
```
