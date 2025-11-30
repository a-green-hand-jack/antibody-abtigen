"""
Epitope-centric antibody-antigen dataset pipeline.

This package implements a pipeline for grouping antibody-antigen structures
by epitope structural similarity (using ESM-2 embeddings) rather than
sequence homology.

Quick Start:
    >>> from antibody_abtigen.epitope_pipeline import create_default_config
    >>> config = create_default_config(data_dir="./data")
    >>> # Run pipeline (implementations in later modules)

Public API:
    - Core data structures and interfaces (from core submodule)
    - Pipeline implementations (from top-level modules)
"""

# Re-export core API
from .core import (
    # Data structures
    ChainMapping,
    CleanedStructure,
    EpitopeResidues,
    EpitopeEmbedding,
    EpitopeGroup,
    AlignedComplex,
    PipelineConfig,

    # Interfaces
    StructureCleaner,
    EpitopeExtractor,
    EpitopeEncoder,
    EmbeddingStore,
    EpitopeGrouper,
    StructureAligner,
    PipelineOrchestrator,

    # Configuration
    load_config_from_yaml,
    create_default_config,
    save_config_to_yaml,

    # Exceptions
    EpitopePipelineError,
    StructureCleaningError,
    EpitopeExtractionError,
    EmbeddingError,
    GroupingError,
    AlignmentError,
    ConfigurationError,
    IndexMappingError,
)

# Implementation modules will be imported here as they're created:
# from .cleaner import GemmiStructureCleaner
# from .extractor import GeometricEpitopeExtractor
# from .encoder import ESM2EpitopeEncoder
# from .storage import HDF5EmbeddingStore
# from .grouper import FAISSEpitopeGrouper
# from .aligner import PyMOLStructureAligner
# from .orchestrator import EpitopePipeline

__version__ = "0.1.0"

__all__ = [
    # Data structures
    'ChainMapping',
    'CleanedStructure',
    'EpitopeResidues',
    'EpitopeEmbedding',
    'EpitopeGroup',
    'AlignedComplex',
    'PipelineConfig',

    # Interfaces
    'StructureCleaner',
    'EpitopeExtractor',
    'EpitopeEncoder',
    'EmbeddingStore',
    'EpitopeGrouper',
    'StructureAligner',
    'PipelineOrchestrator',

    # Configuration
    'load_config_from_yaml',
    'create_default_config',
    'save_config_to_yaml',

    # Exceptions
    'EpitopePipelineError',
    'StructureCleaningError',
    'EpitopeExtractionError',
    'EmbeddingError',
    'GroupingError',
    'AlignmentError',
    'ConfigurationError',
    'IndexMappingError',
]
