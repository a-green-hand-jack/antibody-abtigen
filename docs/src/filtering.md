# filtering.py - Interaction Filtering Module

## Overview

Post-processing filter to validate antibody-antigen interactions by checking atom-atom distances.

## Location

`antibody_abtigen/filtering.py`

## Purpose

After the main pipeline runs, some data points may not have actual antibody-antigen contacts (e.g., if the antigen chain in the PDB is not in contact with the antibody). This module filters out such cases.

## Public API

### Functions

#### `filter_dataset(input_dir, output_dir, distance_threshold=5.0, min_contacts=10, dry_run=False) -> List[FilterResult]`

Filter dataset to keep only data points with sufficient antibody-antigen contacts.

**Parameters:**
- `input_dir`: Directory containing DP_* folders
- `output_dir`: Directory for filtered output
- `distance_threshold`: Max atom-atom distance to count as contact (Å)
- `min_contacts`: Minimum contact pairs required
- `dry_run`: If True, only analyze without copying

**Returns:** List of FilterResult objects

### Classes

#### `FilterResult`

Result of filtering a single data point.

**Fields:**
- `dp_id`: Data point ID
- `passed`: Whether it passed the filter
- `contact_count`: Number of contacts found
- `error`: Error message if failed

## Algorithm

1. Load antibody structure (H + L chains)
2. Load human antigen structure
3. For each atom pair (antibody atom, antigen atom):
   - Calculate Euclidean distance
   - If distance ≤ threshold, count as contact
4. If contact_count ≥ min_contacts, pass filter

## Usage Example

```python
from antibody_abtigen.filtering import filter_dataset

results = filter_dataset(
    input_dir="./output",
    output_dir="./output_filtered",
    distance_threshold=5.0,
    min_contacts=10,
    dry_run=False
)

passed = sum(1 for r in results if r.passed)
print(f"Passed: {passed}/{len(results)}")
```

## CLI Usage

```bash
antibody-abtigen filter-interactions \
  --input ./output \
  --output ./output_filtered \
  --distance-threshold 5.0 \
  --min-contacts 10
```

## Output

- Copies passing DP_* folders to output directory
- Creates `filter_summary.csv` with:
  - `dp_id`: Data point ID
  - `passed`: True/False
  - `contact_count`: Number of contacts
  - `error`: Error message if any
