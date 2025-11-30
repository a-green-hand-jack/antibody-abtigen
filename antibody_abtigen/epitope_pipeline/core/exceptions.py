"""
Custom exceptions for the epitope-centric pipeline.
"""


class EpitopePipelineError(Exception):
    """Base exception for all epitope pipeline errors."""
    pass


class StructureCleaningError(EpitopePipelineError):
    """Raised when structure cleaning fails."""
    pass


class EpitopeExtractionError(EpitopePipelineError):
    """Raised when epitope extraction fails."""
    pass


class EmbeddingError(EpitopePipelineError):
    """Raised when ESM-2 embedding generation fails."""
    pass


class GroupingError(EpitopePipelineError):
    """Raised when similarity grouping fails."""
    pass


class AlignmentError(EpitopePipelineError):
    """Raised when structure alignment fails."""
    pass


class ConfigurationError(EpitopePipelineError):
    """Raised when configuration is invalid."""
    pass


class IndexMappingError(EpitopePipelineError):
    """
    Raised when residue index mapping fails.

    This is critical because we need to map between:
    - ESM-2 0-based indices
    - PDB auth_seq_id
    - PyMOL resi selections
    """
    pass
