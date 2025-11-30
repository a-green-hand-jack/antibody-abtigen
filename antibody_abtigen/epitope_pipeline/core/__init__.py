"""
Core module for the epitope-centric pipeline.

Public API exports:
- Data structures (dataclasses)
- Abstract interfaces
- Configuration management
- Custom exceptions
"""

# Data structures
from .dataclasses import (
    ChainMapping,
    CleanedStructure,
    EpitopeResidues,
    EpitopeEmbedding,
    EpitopeGroup,
    AlignedComplex,
    PipelineConfig
)

# Interfaces
from .interfaces import (
    StructureCleaner,
    EpitopeExtractor,
    EpitopeEncoder,
    EmbeddingStore,
    EpitopeGrouper,
    StructureAligner,
    PipelineOrchestrator
)

# Configuration
from .config import (
    load_config_from_yaml,
    create_default_config,
    save_config_to_yaml
)

# Exceptions
from .exceptions import (
    EpitopePipelineError,
    StructureCleaningError,
    EpitopeExtractionError,
    EmbeddingError,
    GroupingError,
    AlignmentError,
    ConfigurationError,
    IndexMappingError
)

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
