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

# Implementation modules
from .cleaner import (
    GemmiStructureCleaner,
    FilterReason,
    FilterResult,
    save_filter_log,
    generate_filter_summary,
)
from .extractor import GeometricEpitopeExtractor
from .encoder import ESM2EpitopeEncoder, EncoderOutput, ChainEmbeddingResult
from .storage import HDF5EmbeddingStore
from .epitope_log import (
    save_epitope_residues_csv,
    save_epitope_summary_csv,
    save_embedding_stats_csv,
    generate_epitope_report,
)
from .grouper import (
    NumpyEpitopeGrouper,
    GroupMember,
    GroupResult,
    GroupingOutput,
    save_groups_json,
    save_grouping_stats_csv,
    generate_grouping_report,
)
from .aligner import (
    PyMOLStructureAligner,
    AlignmentResult,
    GroupAlignmentOutput,
    save_alignment_summary_csv,
    generate_alignment_report,
)
from .orchestrator import (
    EpitopePipeline,
    PipelineResult,
    PipelineCheckpoint,
    run_pipeline,
)
from .validator import (
    PocketCropper,
    PocketResidues,
    PocketCoordinates,
    StructureValidator,
    StructuralValidationResult,
    ValidatedGroupResult,
    ValidationOutput,
    save_validated_groups_json,
    save_validation_report_csv,
    generate_validation_report,
    load_validated_groups_json,
)

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

    # Implementations
    'GemmiStructureCleaner',
    'GeometricEpitopeExtractor',
    'ESM2EpitopeEncoder',
    'EncoderOutput',
    'ChainEmbeddingResult',
    'HDF5EmbeddingStore',

    # Filtering
    'FilterReason',
    'FilterResult',
    'save_filter_log',
    'generate_filter_summary',

    # Epitope logging
    'save_epitope_residues_csv',
    'save_epitope_summary_csv',
    'save_embedding_stats_csv',
    'generate_epitope_report',

    # Grouping
    'NumpyEpitopeGrouper',
    'GroupMember',
    'GroupResult',
    'GroupingOutput',
    'save_groups_json',
    'save_grouping_stats_csv',
    'generate_grouping_report',

    # Alignment
    'PyMOLStructureAligner',
    'AlignmentResult',
    'GroupAlignmentOutput',
    'save_alignment_summary_csv',
    'generate_alignment_report',

    # Orchestrator
    'EpitopePipeline',
    'PipelineResult',
    'PipelineCheckpoint',
    'run_pipeline',

    # Structure Validation
    'PocketCropper',
    'PocketResidues',
    'PocketCoordinates',
    'StructureValidator',
    'StructuralValidationResult',
    'ValidatedGroupResult',
    'ValidationOutput',
    'save_validated_groups_json',
    'save_validation_report_csv',
    'generate_validation_report',
    'load_validated_groups_json',
]
