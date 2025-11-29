# cli.py - Command-Line Interface Module

## Overview

Click-based CLI providing three subcommands for building datasets, converting formats, and filtering results.

## Location

`antibody_abtigen/cli.py`

## Entry Point

```bash
# Via installed package
antibody-abtigen [COMMAND] [OPTIONS]

# Via uv run
uv run antibody-abtigen [COMMAND] [OPTIONS]

# Via run.py (defaults to 'build')
uv run python run.py [OPTIONS]
```

## Commands

### `build`

Build cross-species antibody-antigen structure dataset.

```bash
antibody-abtigen build [OPTIONS]
```

**Options:**
| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--output` | `-o` | `./output` | Output directory |
| `--data-dir` | `-d` | `./data` | Cache directory |
| `--limit` | `-n` | All | Max entries to process |
| `--resolution` | `-r` | 2.5 | Max resolution (Å) |
| `--identity` | `-i` | 50.0 | Min sequence identity (%) |
| `--dry-run` | | False | Analysis only |
| `--no-pymol` | | False | Use Biopython alignment |
| `--force-download` | | False | Re-download SAbDab |

**Examples:**
```bash
# Quick test with dry run
antibody-abtigen build --limit 10 --dry-run

# Full run with custom thresholds
antibody-abtigen build --limit 100 --resolution 3.0 --identity 60

# Use Biopython alignment
antibody-abtigen build --limit 5 --no-pymol
```

### `to-yaml`

Convert CIF files to Boltz-1 YAML format.

```bash
antibody-abtigen to-yaml --input DIR --output DIR [OPTIONS]
```

**Options:**
| Option | Short | Required | Description |
|--------|-------|----------|-------------|
| `--input` | `-i` | Yes | Input directory with DP_* folders |
| `--output` | `-o` | Yes | Output directory for YAML files |
| `--no-bonds` | | No | Exclude disulfide bonds |

**Example:**
```bash
antibody-abtigen to-yaml --input ./output --output ./yamls
```

### `filter-interactions`

Filter data points by antibody-antigen contacts.

```bash
antibody-abtigen filter-interactions --input DIR --output DIR [OPTIONS]
```

**Options:**
| Option | Default | Description |
|--------|---------|-------------|
| `--input` | Required | Input directory with DP_* folders |
| `--output` | Required | Output directory |
| `--distance-threshold` | 5.0 | Max contact distance (Å) |
| `--min-contacts` | 10 | Min required contacts |
| `--dry-run` | False | Analysis only |

**Example:**
```bash
antibody-abtigen filter-interactions \
  --input ./output \
  --output ./filtered \
  --distance-threshold 5.0 \
  --min-contacts 10
```

## Usage Workflow

Typical workflow:

```bash
# 1. Build dataset (dry run first)
antibody-abtigen build --limit 10 --dry-run

# 2. Build actual dataset
antibody-abtigen build --limit 100

# 3. Filter for contacts
antibody-abtigen filter-interactions \
  --input ./output \
  --output ./output_filtered

# 4. Convert to YAML for Boltz
antibody-abtigen to-yaml \
  --input ./output_filtered \
  --output ./yamls
```

## Exit Codes

- `0`: Success
- `1`: Error (see stderr for details)
