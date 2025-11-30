"""
Configuration management for the epitope-centric pipeline.

Supports loading from YAML files and environment variables.
"""

import os
from pathlib import Path
from typing import Dict, Optional
import yaml

from .dataclasses import PipelineConfig
from .exceptions import ConfigurationError


def load_config_from_yaml(yaml_path: Path) -> PipelineConfig:
    """
    Load pipeline configuration from a YAML file.

    Args:
        yaml_path: Path to YAML configuration file

    Returns:
        PipelineConfig object

    Raises:
        ConfigurationError: If YAML is invalid or config validation fails

    Example YAML:
        data:
          raw_cif_dir: data/raw_cif
          cleaned_cif_dir: data/epitope_pipeline/cleaned
          aligned_cif_dir: data/epitope_pipeline/aligned_cif
          embeddings_dir: data/epitope_pipeline/embeddings
          meta_dir: data/epitope_pipeline/meta

        epitope_extraction:
          contact_distance_threshold: 5.0

        embedding:
          model_name: esm2_t36_3B_UR50D
          batch_size: 8
          device: cuda
          use_fp16: true

        grouping:
          similarity_threshold: 0.85
          top_k_neighbors: 100

        alignment:
          method: pymol
          pymol_env: /path/to/pymol-env
          use_super: true
    """
    if not yaml_path.exists():
        raise ConfigurationError(f"Config file not found: {yaml_path}")

    try:
        with open(yaml_path, 'r') as f:
            raw_config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigurationError(f"Failed to parse YAML: {e}")

    # Flatten nested structure
    config_dict = _flatten_yaml_config(raw_config)

    # Create config object
    try:
        config = PipelineConfig.from_dict(config_dict)
    except TypeError as e:
        raise ConfigurationError(f"Invalid config parameters: {e}")

    # Validate
    errors = config.validate()
    if errors:
        raise ConfigurationError(f"Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors))

    return config


def _flatten_yaml_config(raw_config: Dict) -> Dict:
    """
    Flatten nested YAML structure to match PipelineConfig fields.

    Converts:
        {
          "data": {"raw_cif_dir": "..."},
          "embedding": {"device": "cuda"}
        }
    To:
        {
          "raw_cif_dir": "...",
          "device": "cuda"
        }
    """
    flattened = {}

    # Data directories
    if 'data' in raw_config:
        data_section = raw_config['data']
        # Ensure data_dir is set (required base path)
        if 'data_dir' not in data_section:
            # Infer from other paths if possible
            if 'raw_cif_dir' in data_section:
                # Try to infer base data_dir
                raw_path = Path(data_section['raw_cif_dir'])
                if 'raw_cif' in str(raw_path):
                    data_section['data_dir'] = str(raw_path.parent)
                else:
                    data_section['data_dir'] = './data/epitope_pipeline'
            else:
                data_section['data_dir'] = './data/epitope_pipeline'

        for key, value in data_section.items():
            flattened[key] = value

    # Epitope extraction
    if 'epitope_extraction' in raw_config:
        for key, value in raw_config['epitope_extraction'].items():
            flattened[key] = value

    # Embedding
    if 'embedding' in raw_config:
        embed_section = raw_config['embedding']
        flattened['esm_model_name'] = embed_section.get('model_name', 'esm2_t36_3B_UR50D')
        flattened['embedding_batch_size'] = embed_section.get('batch_size', 8)
        flattened['device'] = embed_section.get('device', 'cuda')
        flattened['use_fp16'] = embed_section.get('use_fp16', True)
        if 'cache_dir' in embed_section:
            flattened['esm_cache_dir'] = embed_section['cache_dir']

    # Grouping
    if 'grouping' in raw_config:
        group_section = raw_config['grouping']
        flattened['similarity_metric'] = group_section.get('similarity_metric', 'cosine')
        flattened['similarity_threshold'] = group_section.get('similarity_threshold', 0.85)
        flattened['top_k_neighbors'] = group_section.get('top_k_neighbors', 100)

    # Alignment
    if 'alignment' in raw_config:
        align_section = raw_config['alignment']
        flattened['alignment_method'] = align_section.get('method', 'pymol')
        flattened['use_pymol_super'] = align_section.get('use_super', True)
        if 'pymol_env' in align_section:
            flattened['pymol_env_path'] = align_section['pymol_env']

    # Performance
    if 'performance' in raw_config:
        perf_section = raw_config['performance']
        flattened['parallel_downloads'] = perf_section.get('parallel_downloads', 16)
        flattened['parallel_alignments'] = perf_section.get('parallel_alignments', 12)
        flattened['gpu_memory_fraction'] = perf_section.get('gpu_memory_fraction', 0.9)

    # Clustering
    if 'clustering' in raw_config:
        clust_section = raw_config['clustering']
        flattened['antigen_clustering_identity'] = clust_section.get('antigen_identity_threshold', 0.90)

    return flattened


def create_default_config(
    data_dir: Path,
    device: str = "cuda",
    use_fp16: bool = True,
    similarity_threshold: float = 0.85,
    pymol_env_path: Optional[Path] = None
) -> PipelineConfig:
    """
    Create a default configuration with sensible defaults.

    Args:
        data_dir: Base data directory
        device: 'cuda' or 'cpu'
        use_fp16: Use FP16 mixed precision for ESM-2
        similarity_threshold: Cosine similarity threshold for grouping
        pymol_env_path: Path to PyMOL conda environment (optional)

    Returns:
        PipelineConfig with default settings
    """
    data_dir = Path(data_dir)

    config = PipelineConfig(
        data_dir=data_dir,
        raw_cif_dir=data_dir / "raw_cif",
        cleaned_cif_dir=data_dir / "epitope_pipeline" / "cleaned",
        aligned_cif_dir=data_dir / "epitope_pipeline" / "aligned_cif",
        embeddings_dir=data_dir / "epitope_pipeline" / "embeddings",
        meta_dir=data_dir / "epitope_pipeline" / "meta",
        contact_distance_threshold=5.0,
        esm_model_name="esm2_t36_3B_UR50D",
        embedding_batch_size=8,
        device=device,
        use_fp16=use_fp16,
        similarity_metric="cosine",
        similarity_threshold=similarity_threshold,
        top_k_neighbors=100,
        alignment_method="pymol",
        pymol_env_path=pymol_env_path,
        use_pymol_super=True,
        parallel_downloads=16,
        parallel_alignments=12,
        gpu_memory_fraction=0.9,
        antigen_clustering_identity=0.90
    )

    # Validate
    errors = config.validate()
    if errors:
        raise ConfigurationError(f"Default config validation failed:\n" + "\n".join(f"  - {e}" for e in errors))

    return config


def save_config_to_yaml(config: PipelineConfig, output_path: Path) -> None:
    """
    Save configuration to a YAML file.

    Args:
        config: PipelineConfig to save
        output_path: Where to save YAML file
    """
    # Convert to nested structure
    yaml_dict = {
        'data': {
            'data_dir': str(config.data_dir),
            'raw_cif_dir': str(config.raw_cif_dir),
            'cleaned_cif_dir': str(config.cleaned_cif_dir),
            'aligned_cif_dir': str(config.aligned_cif_dir),
            'embeddings_dir': str(config.embeddings_dir),
            'meta_dir': str(config.meta_dir),
        },
        'epitope_extraction': {
            'contact_distance_threshold': config.contact_distance_threshold,
        },
        'embedding': {
            'model_name': config.esm_model_name,
            'batch_size': config.embedding_batch_size,
            'device': config.device,
            'use_fp16': config.use_fp16,
        },
        'grouping': {
            'similarity_metric': config.similarity_metric,
            'similarity_threshold': config.similarity_threshold,
            'top_k_neighbors': config.top_k_neighbors,
        },
        'alignment': {
            'method': config.alignment_method,
            'use_super': config.use_pymol_super,
        },
        'performance': {
            'parallel_downloads': config.parallel_downloads,
            'parallel_alignments': config.parallel_alignments,
            'gpu_memory_fraction': config.gpu_memory_fraction,
        },
        'clustering': {
            'antigen_identity_threshold': config.antigen_clustering_identity,
        }
    }

    # Add optional fields
    if config.esm_cache_dir:
        yaml_dict['embedding']['cache_dir'] = str(config.esm_cache_dir)

    if config.pymol_env_path:
        yaml_dict['alignment']['pymol_env'] = str(config.pymol_env_path)

    # Write to file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False, indent=2)
