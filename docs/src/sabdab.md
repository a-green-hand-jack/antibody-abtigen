# sabdab.py - SAbDab Data Module

## Overview

Handles downloading and parsing data from the Structural Antibody Database (SAbDab) hosted at Oxford University.

## Location

`antibody_abtigen/sabdab.py`

## Public API

### Functions

#### `download_sabdab_summary(output_dir: str, force_download: bool = False) -> str`

Download the SAbDab summary TSV file.

**Parameters:**
- `output_dir`: Directory to save the summary file
- `force_download`: If True, re-download even if file exists

**Returns:** Path to the downloaded summary file

**Example:**
```python
from antibody_abtigen.sabdab import download_sabdab_summary

summary_path = download_sabdab_summary("./data")
# Returns: "./data/sabdab_summary.tsv"
```

#### `parse_sabdab_summary(summary_path: str) -> pd.DataFrame`

Parse the SAbDab summary TSV file into a DataFrame.

**Parameters:**
- `summary_path`: Path to the summary file

**Returns:** DataFrame with parsed data

#### `filter_human_antigen_complexes(df, resolution_threshold=2.5, require_paired=True) -> pd.DataFrame`

Filter entries for human antigen protein complexes.

**Filtering steps:**
1. Human antigens only (`antigen_species` contains "homo sapiens")
2. Protein antigens only (`antigen_type` contains "protein")
3. Resolution ≤ threshold (default 2.5Å)
4. Paired H/L chains present (if `require_paired=True`)
5. Antigen chain information available

**Parameters:**
- `df`: SAbDab summary DataFrame
- `resolution_threshold`: Maximum resolution in Angstroms
- `require_paired`: Require both heavy and light chains

**Returns:** Filtered DataFrame

#### `get_unique_complexes(df: pd.DataFrame) -> pd.DataFrame`

Remove duplicate complexes based on PDB ID and antigen chain.

#### `extract_chain_info(row: pd.Series) -> dict`

Extract chain information from a single SAbDab entry.

## Constants

- `SABDAB_SUMMARY_URL`: URL to SAbDab summary endpoint

## Data Source

- **URL**: https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/
- **Format**: Tab-separated values (TSV)
- **Updates**: Database is regularly updated by Oxford

## Usage Example

```python
from antibody_abtigen.sabdab import (
    download_sabdab_summary,
    parse_sabdab_summary,
    filter_human_antigen_complexes,
    get_unique_complexes
)

# Download and parse
summary_path = download_sabdab_summary("./data")
df = parse_sabdab_summary(summary_path)

# Filter
filtered = filter_human_antigen_complexes(df, resolution_threshold=2.5)
unique = get_unique_complexes(filtered)

print(f"Found {len(unique)} unique human antigen complexes")
```
