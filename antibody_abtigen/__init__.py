"""
Antibody-Antigen Cross-Species Structure Dataset Builder

A pipeline for building cross-species antibody-antigen structure datasets
from SAbDab and PDB data.
"""

__version__ = "0.1.0"

from .pipeline import CrossSpeciesDatasetPipeline
from .sabdab import download_sabdab_summary, filter_human_antigen_complexes
from .mapping import (
    get_uniprot_from_pdb_chain,
    find_mouse_ortholog,
    calculate_sequence_identity,
)
from .structure import (
    download_structure,
    parse_structure,
    align_structures_biopython,
    align_structures_pymol,
)

__all__ = [
    "CrossSpeciesDatasetPipeline",
    "download_sabdab_summary",
    "filter_human_antigen_complexes",
    "get_uniprot_from_pdb_chain",
    "find_mouse_ortholog",
    "calculate_sequence_identity",
    "download_structure",
    "parse_structure",
    "align_structures_biopython",
    "align_structures_pymol",
]
