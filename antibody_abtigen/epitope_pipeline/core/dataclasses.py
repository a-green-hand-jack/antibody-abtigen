"""
Core data structures for the epitope-centric pipeline.

These immutable dataclasses represent the state at each stage of the pipeline:
1. CleanedStructure - After structure cleaning (Phase 1)
2. EpitopeResidues - After epitope extraction (Phase 2)
3. EpitopeEmbedding - After ESM-2 encoding (Phase 3)
4. EpitopeGroup - After similarity grouping (Phase 4)
5. AlignedComplex - After PyMOL alignment (Phase 5)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class ChainMapping:
    """
    Maps original PDB chain IDs to standardized chain IDs.

    This is crucial for index mapping between different systems:
    - ESM-2 uses 0-based sequence indices
    - PyMOL uses PDB auth_seq_id
    - We need to track both to convert between them
    """
    pdb_id: str
    original_chain_id: str
    standardized_chain_id: str  # A, B for antigen; H, L for antibody
    chain_type: str  # 'antigen' | 'antibody_heavy' | 'antibody_light'
    sequence: str  # Full amino acid sequence
    auth_seq_id_map: Dict[int, int] = field(default_factory=dict)  # 0-based -> PDB auth_seq_id


@dataclass(frozen=True)
class CleanedStructure:
    """
    Result of structure cleaning (Phase 1).

    Contains:
    - Cleaned CIF file path
    - Chain ID mappings (original -> standardized)
    - Antigen cluster assignment (for filtering same-antigen pairs)
    """
    pdb_id: str
    file_path: Path
    chain_mappings: List[ChainMapping]
    antigen_cluster_id: Optional[str] = None
    metadata: Dict = field(default_factory=dict)

    def get_chain_mapping(self, standardized_id: str) -> Optional[ChainMapping]:
        """Get chain mapping by standardized chain ID (e.g., 'A', 'H')."""
        for mapping in self.chain_mappings:
            if mapping.standardized_chain_id == standardized_id:
                return mapping
        return None

    def get_antigen_chains(self) -> List[ChainMapping]:
        """Get all antigen chain mappings."""
        return [m for m in self.chain_mappings if m.chain_type == 'antigen']

    def get_antibody_chains(self) -> List[ChainMapping]:
        """Get all antibody chain mappings."""
        return [m for m in self.chain_mappings
                if m.chain_type in ('antibody_heavy', 'antibody_light')]


@dataclass(frozen=True)
class EpitopeResidues:
    """
    Defines epitope residues in a structure (Phase 2).

    Epitope is defined as antigen residues within distance_threshold
    of any antibody atom. Supports multi-chain epitopes (e.g., when
    antibody binds at the interface of a trimer).

    Residue numbers use PDB auth_seq_id (not 0-based indices).
    """
    epitope_id: str  # Unique identifier (e.g., "7k8t_epi")
    pdb_id: str
    antigen_chains: Dict[str, List[int]]  # {"A": [100,101,102], "B": [45,46]}
    antibody_chains: List[str]  # ["H", "L"]
    distance_threshold: float = 5.0
    num_contacts: int = 0  # Total number of atom-atom contacts

    def get_all_residues(self) -> List[Tuple[str, int]]:
        """Get flattened list of all epitope residues as (chain_id, res_num)."""
        residues = []
        for chain_id, res_nums in self.antigen_chains.items():
            for res_num in res_nums:
                residues.append((chain_id, res_num))
        return sorted(residues)

    def total_residue_count(self) -> int:
        """Get total number of epitope residues across all chains."""
        return sum(len(res_list) for res_list in self.antigen_chains.values())


@dataclass(frozen=True)
class EpitopeEmbedding:
    """
    ESM-2 embedding for an epitope (Phase 3).

    Uses full-context embedding strategy:
    1. Process entire antigen sequence through ESM-2
    2. Extract embeddings for epitope residue positions
    3. Mean pooling to get single 2560-dim vector
    4. L2 normalization for cosine similarity
    """
    epitope_id: str
    pdb_id: str
    embedding: np.ndarray  # (2560,) float32, L2-normalized
    full_sequence: str  # Complete antigen sequence used for context
    epitope_residue_indices: List[int]  # 0-based indices into full_sequence
    embedding_method: str = "esm2_t36_3B_UR50D"
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        """Validate embedding shape and normalization."""
        if self.embedding.shape != (2560,):
            raise ValueError(f"Embedding must be shape (2560,), got {self.embedding.shape}")
        # Check L2 normalization (allow small floating point error)
        norm = np.linalg.norm(self.embedding)
        if not np.isclose(norm, 1.0, atol=1e-5):
            raise ValueError(f"Embedding must be L2-normalized, got norm={norm:.6f}")


@dataclass(frozen=True)
class EpitopeGroup:
    """
    A group of structurally similar epitopes (Phase 4).

    Epitopes are grouped by:
    1. High cosine similarity (>= threshold, e.g., 0.85)
    2. Different antigen clusters (avoid grouping same antigen)
    3. Connected components in similarity graph (NetworkX)
    """
    group_id: str  # e.g., "group_001"
    reference_epitope_id: str  # Anchor for alignment
    member_epitope_ids: List[str]  # All members including reference
    avg_similarity: float  # Average pairwise cosine similarity
    metadata: Dict = field(default_factory=dict)

    def size(self) -> int:
        """Get number of members in this group."""
        return len(self.member_epitope_ids)

    def get_mobile_ids(self) -> List[str]:
        """Get epitope IDs that need to be aligned (exclude reference)."""
        return [eid for eid in self.member_epitope_ids if eid != self.reference_epitope_id]


@dataclass(frozen=True)
class AlignedComplex:
    """
    Result of PyMOL epitope-based alignment (Phase 5).

    Contains the aligned full complex (antigen + antibody) in CIF format.
    Alignment is based on epitope residues, but transformation is applied
    to the entire complex to preserve relative antibody position.
    """
    epitope_id: str
    pdb_id: str
    group_id: str
    aligned_file_path: Path
    rmsd_to_reference: float  # RMSD of epitope CA atoms after alignment
    alignment_method: str  # 'pymol_super' | 'biopython'
    is_reference: bool = False
    metadata: Dict = field(default_factory=dict)


@dataclass
class PipelineConfig:
    """
    Configuration for the epitope-centric pipeline.

    This is mutable to allow runtime adjustments, but should generally
    be loaded from a YAML file and kept constant during a run.
    """
    # Directory structure
    data_dir: Path
    raw_cif_dir: Path
    cleaned_cif_dir: Path
    aligned_cif_dir: Path
    embeddings_dir: Path
    meta_dir: Path

    # Epitope extraction parameters
    contact_distance_threshold: float = 5.0

    # Embedding parameters
    esm_model_name: str = "esm2_t36_3B_UR50D"
    embedding_batch_size: int = 8
    device: str = "cuda"
    use_fp16: bool = True
    esm_cache_dir: Optional[Path] = None

    # Grouping parameters
    similarity_metric: str = "cosine"
    similarity_threshold: float = 0.85
    top_k_neighbors: int = 100

    # Alignment parameters
    alignment_method: str = "pymol"  # 'pymol' | 'biopython'
    pymol_env_path: Optional[Path] = None
    use_pymol_super: bool = True

    # Performance parameters
    parallel_downloads: int = 16
    parallel_alignments: int = 12
    gpu_memory_fraction: float = 0.9

    # Clustering parameters (for antigen deduplication)
    antigen_clustering_identity: float = 0.90

    def validate(self) -> List[str]:
        """
        Validate configuration and return list of errors.

        Returns:
            Empty list if valid, otherwise list of error messages
        """
        errors = []

        # Check directories exist (or can be created)
        for dir_attr in ['data_dir', 'raw_cif_dir', 'cleaned_cif_dir',
                         'aligned_cif_dir', 'embeddings_dir', 'meta_dir']:
            dir_path = getattr(self, dir_attr)
            if not isinstance(dir_path, Path):
                errors.append(f"{dir_attr} must be a Path object, got {type(dir_path)}")

        # Check thresholds are in valid ranges
        if not 0.0 < self.contact_distance_threshold <= 10.0:
            errors.append(f"contact_distance_threshold must be in (0, 10], got {self.contact_distance_threshold}")

        if not 0.0 <= self.similarity_threshold <= 1.0:
            errors.append(f"similarity_threshold must be in [0, 1], got {self.similarity_threshold}")

        if not 0.0 <= self.antigen_clustering_identity <= 1.0:
            errors.append(f"antigen_clustering_identity must be in [0, 1], got {self.antigen_clustering_identity}")

        # Check device
        if self.device not in ['cuda', 'cpu']:
            errors.append(f"device must be 'cuda' or 'cpu', got {self.device}")

        # Check alignment method
        if self.alignment_method not in ['pymol', 'biopython']:
            errors.append(f"alignment_method must be 'pymol' or 'biopython', got {self.alignment_method}")

        return errors

    @classmethod
    def from_dict(cls, config_dict: Dict) -> 'PipelineConfig':
        """
        Create config from dictionary (e.g., loaded from YAML).

        Converts string paths to Path objects.
        """
        # Convert string paths to Path objects
        path_fields = ['data_dir', 'raw_cif_dir', 'cleaned_cif_dir',
                      'aligned_cif_dir', 'embeddings_dir', 'meta_dir',
                      'esm_cache_dir', 'pymol_env_path']

        for field_name in path_fields:
            if field_name in config_dict and config_dict[field_name] is not None:
                config_dict[field_name] = Path(config_dict[field_name])

        return cls(**config_dict)
