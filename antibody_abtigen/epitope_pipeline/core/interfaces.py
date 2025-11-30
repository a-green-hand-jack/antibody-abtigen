"""
Abstract base classes defining interfaces for pipeline components.

This enables:
1. Dependency injection for testing
2. Easy swapping of implementations (e.g., different embedding models)
3. Clear separation of concerns
"""

from abc import ABC, abstractmethod
from typing import List, Iterator, Tuple, Dict, Optional
from pathlib import Path

from .dataclasses import (
    CleanedStructure,
    EpitopeResidues,
    EpitopeEmbedding,
    EpitopeGroup,
    AlignedComplex,
    PipelineConfig
)


class StructureCleaner(ABC):
    """
    Interface for structure cleaning operations.

    Implementations:
    - GemmiStructureCleaner: Uses gemmi library for mmCIF processing
    """

    @abstractmethod
    def clean_structure(self, cif_path: Path) -> CleanedStructure:
        """
        Clean and standardize a single structure.

        Args:
            cif_path: Path to raw mmCIF file

        Returns:
            CleanedStructure with standardized chain IDs and index mappings

        Raises:
            StructureCleaningError: If cleaning fails
        """
        pass

    @abstractmethod
    def batch_clean(self, cif_paths: List[Path]) -> Iterator[CleanedStructure]:
        """
        Clean multiple structures (generator for memory efficiency).

        Args:
            cif_paths: List of paths to raw mmCIF files

        Yields:
            CleanedStructure objects (may skip failed structures)
        """
        pass


class EpitopeExtractor(ABC):
    """
    Interface for epitope extraction.

    Implementations:
    - GeometricEpitopeExtractor: Distance-based contact detection
    """

    @abstractmethod
    def extract_epitope(
        self,
        structure: CleanedStructure,
        distance_threshold: float
    ) -> EpitopeResidues:
        """
        Extract epitope residues from a structure.

        Args:
            structure: Cleaned structure
            distance_threshold: Contact distance cutoff in Angstroms

        Returns:
            EpitopeResidues defining the epitope

        Raises:
            EpitopeExtractionError: If extraction fails
        """
        pass


class EpitopeEncoder(ABC):
    """
    Interface for epitope embedding.

    Implementations:
    - ESM2EpitopeEncoder: Uses ESM-2 3B model with full-context strategy
    - ProtTransEncoder: Alternative using ProtTrans (for future extension)
    """

    @abstractmethod
    def encode(
        self,
        epitope: EpitopeResidues,
        structure: CleanedStructure
    ) -> EpitopeEmbedding:
        """
        Generate embedding for a single epitope.

        Args:
            epitope: Epitope residue definition
            structure: Cleaned structure containing the antigen

        Returns:
            EpitopeEmbedding with normalized vector

        Raises:
            EmbeddingError: If encoding fails
        """
        pass

    @abstractmethod
    def batch_encode(
        self,
        epitopes: List[Tuple[EpitopeResidues, CleanedStructure]]
    ) -> List[EpitopeEmbedding]:
        """
        Batch encode for GPU efficiency.

        Args:
            epitopes: List of (epitope, structure) pairs

        Returns:
            List of EpitopeEmbedding objects

        Raises:
            EmbeddingError: If batch encoding fails
        """
        pass


class EmbeddingStore(ABC):
    """
    Interface for storing and retrieving embeddings.

    Implementations:
    - HDF5EmbeddingStore: Uses HDF5 with gzip compression
    """

    @abstractmethod
    def save(self, embeddings: List[EpitopeEmbedding], output_path: Path) -> None:
        """
        Save embeddings to storage.

        Args:
            embeddings: List of embeddings to save
            output_path: Where to save (e.g., HDF5 file path)
        """
        pass

    @abstractmethod
    def load(self, input_path: Path) -> List[EpitopeEmbedding]:
        """
        Load embeddings from storage.

        Args:
            input_path: Where to load from

        Returns:
            List of loaded embeddings
        """
        pass

    @abstractmethod
    def get_by_id(self, input_path: Path, epitope_id: str) -> Optional[EpitopeEmbedding]:
        """
        Load a single embedding by ID.

        Args:
            input_path: Storage location
            epitope_id: ID of epitope to retrieve

        Returns:
            EpitopeEmbedding if found, None otherwise
        """
        pass


class EpitopeGrouper(ABC):
    """
    Interface for grouping epitopes by similarity.

    Implementations:
    - FAISSEpitopeGrouper: Uses FAISS + NetworkX for fast GPU-accelerated grouping
    """

    @abstractmethod
    def build_index(self, embeddings: List[EpitopeEmbedding]) -> None:
        """
        Build search index from embeddings.

        Args:
            embeddings: All epitope embeddings to index

        Raises:
            GroupingError: If index building fails
        """
        pass

    @abstractmethod
    def find_groups(
        self,
        similarity_threshold: float,
        antigen_clusters: Optional[Dict[str, str]] = None
    ) -> List[EpitopeGroup]:
        """
        Find connected components of similar epitopes.

        Args:
            similarity_threshold: Minimum similarity to create edge
            antigen_clusters: Mapping epitope_id -> cluster_id (for filtering)

        Returns:
            List of EpitopeGroup objects

        Raises:
            GroupingError: If grouping fails
        """
        pass


class StructureAligner(ABC):
    """
    Interface for structure alignment.

    Implementations:
    - PyMOLStructureAligner: Uses PyMOL super command (preferred)
    - BiopythonStructureAligner: Fallback using Biopython Superimposer
    """

    @abstractmethod
    def align_group(
        self,
        group: EpitopeGroup,
        epitopes: Dict[str, EpitopeResidues],
        structures: Dict[str, CleanedStructure],
        output_dir: Path
    ) -> List[AlignedComplex]:
        """
        Align all structures in a group to the reference.

        Args:
            group: Group definition with reference and members
            epitopes: Mapping epitope_id -> EpitopeResidues
            structures: Mapping pdb_id -> CleanedStructure
            output_dir: Where to save aligned CIF files

        Returns:
            List of AlignedComplex objects (including reference)

        Raises:
            AlignmentError: If alignment fails
        """
        pass


class PipelineOrchestrator(ABC):
    """
    Interface for high-level pipeline orchestration.

    Implementations:
    - EpitopePipeline: Main orchestrator with checkpointing
    """

    @abstractmethod
    def run(
        self,
        input_cif_paths: List[Path],
        output_json: Path,
        checkpoint_every: int = 100
    ) -> Dict:
        """
        Run the full 5-stage pipeline.

        Args:
            input_cif_paths: List of raw CIF files to process
            output_json: Where to save final dataset JSON
            checkpoint_every: Save checkpoint after N structures

        Returns:
            Summary dictionary with statistics

        Raises:
            EpitopePipelineError: If pipeline fails
        """
        pass
