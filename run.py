#!/usr/bin/env python3
"""
Cross-Species Antibody-Antigen Structure Dataset Builder

This is a convenience script for running the pipeline directly.
You can also use the installed CLI: `antibody-abtigen --help`

Usage:
    # Full run with all data
    uv run python run.py --output ./output

    # Demo run with limited entries
    uv run python run.py --output ./output --limit 10

    # Dry run (analysis only, no structure downloads)
    uv run python run.py --output ./output --limit 50 --dry-run

    # Custom thresholds
    uv run python run.py --output ./output --resolution 3.0 --identity 60
"""

import sys

from antibody_abtigen.cli import main

if __name__ == '__main__':
    # Convenience: allow `python run.py --limit 10` to default to `build`
    known_commands = {'build', 'to-yaml', 'filter-interactions', '--help', '--version'}
    if len(sys.argv) == 1:
        sys.argv.insert(1, 'build')
    elif sys.argv[1] not in known_commands:
        sys.argv.insert(1, 'build')

    main()
