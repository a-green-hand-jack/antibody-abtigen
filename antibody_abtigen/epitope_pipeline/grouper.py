"""
Epitope grouping module for the epitope-centric pipeline.

This module implements Phase 4: Similarity-based Grouping.

Key features:
1. Compute pairwise cosine similarity between epitope embeddings
2. Build similarity graph with NetworkX
3. Find connected components as epitope groups
4. Select reference (most connected node) for each group
5. Output groups.json and similarity matrix

For L2-normalized embeddings, cosine similarity = dot product.
"""

import json
import csv
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, field, asdict

import numpy as np
import networkx as nx
import h5py

from .core import (
    EpitopeGrouper,
    EpitopeEmbedding,
    EpitopeGroup,
    GroupingError,
)
from .encoder import EncoderOutput
from .storage import HDF5EmbeddingStore

logger = logging.getLogger(__name__)


@dataclass
class GroupMember:
    """A member of an epitope group."""
    epitope_id: str
    pdb_id: str
    is_reference: bool
    antigen_chains: List[str]
    epitope_residues: Dict[str, List[int]]
    total_residues: int
    similarity_to_ref: Optional[float] = None  # None for reference


@dataclass
class GroupResult:
    """Result of grouping for a single group."""
    group_id: str
    reference_epitope_id: str
    member_count: int
    avg_similarity: float
    min_similarity: float
    max_similarity: float
    members: List[GroupMember]


@dataclass
class GroupingOutput:
    """Complete output from the grouper."""
    groups: List[GroupResult]
    total_epitopes: int
    grouped_epitopes: int
    singleton_count: int  # Epitopes with no similar matches
    similarity_threshold: float

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'metadata': {
                'total_epitopes': self.total_epitopes,
                'grouped_epitopes': self.grouped_epitopes,
                'singleton_count': self.singleton_count,
                'num_groups': len(self.groups),
                'similarity_threshold': self.similarity_threshold,
            },
            'groups': [asdict(g) for g in self.groups]
        }


class NumpyEpitopeGrouper(EpitopeGrouper):
    """
    NumPy + NetworkX based epitope grouper.

    Uses pure NumPy for similarity computation (efficient for L2-normalized vectors)
    and NetworkX for graph-based clustering.

    For L2-normalized embeddings: cosine_similarity = dot_product

    Usage:
        >>> grouper = NumpyEpitopeGrouper(similarity_threshold=0.85)
        >>> grouper.build_index(embeddings)
        >>> groups = grouper.find_groups()
    """

    def __init__(
        self,
        similarity_threshold: float = 0.85,
        min_group_size: int = 2,
        exclude_same_pdb: bool = True
    ):
        """
        Initialize grouper.

        Args:
            similarity_threshold: Minimum cosine similarity to create edge (default: 0.85)
            min_group_size: Minimum members in a group (default: 2)
            exclude_same_pdb: Exclude pairs from same PDB structure (default: True)
        """
        self.similarity_threshold = similarity_threshold
        self.min_group_size = min_group_size
        self.exclude_same_pdb = exclude_same_pdb

        self._embeddings: Optional[np.ndarray] = None
        self._epitope_ids: Optional[List[str]] = None
        self._pdb_ids: Optional[List[str]] = None
        self._epitope_metadata: Optional[Dict[str, Dict]] = None
        self._similarity_matrix: Optional[np.ndarray] = None

    def build_index(
        self,
        embeddings: List[EpitopeEmbedding],
        epitope_metadata: Optional[Dict[str, Dict]] = None
    ) -> None:
        """
        Build index from embeddings.

        Args:
            embeddings: List of EpitopeEmbedding objects
            epitope_metadata: Optional dict mapping epitope_id -> metadata
                             (antigen_chains, epitope_residues, etc.)

        Raises:
            GroupingError: If index building fails
        """
        if not embeddings:
            raise GroupingError("No embeddings provided")

        try:
            # Extract embedding vectors
            self._epitope_ids = [e.epitope_id for e in embeddings]
            self._pdb_ids = [e.pdb_id for e in embeddings]

            # Stack into matrix (N, 2560)
            self._embeddings = np.stack([e.embedding for e in embeddings])

            # Store metadata
            self._epitope_metadata = epitope_metadata or {}

            logger.info(f"Built index with {len(embeddings)} embeddings")

        except Exception as e:
            raise GroupingError(f"Failed to build index: {e}") from e

    def build_index_from_encoder_outputs(
        self,
        outputs: List[EncoderOutput],
        use_epitope_embedding: bool = True
    ) -> None:
        """
        Build index from EncoderOutput objects.

        Args:
            outputs: List of EncoderOutput objects
            use_epitope_embedding: Use epitope embedding (True) or full embedding (False)
        """
        if not outputs:
            raise GroupingError("No encoder outputs provided")

        try:
            self._epitope_ids = [o.epitope_id for o in outputs]
            self._pdb_ids = [o.pdb_id for o in outputs]

            # Choose embedding type
            if use_epitope_embedding:
                self._embeddings = np.stack([o.epitope_embedding for o in outputs])
            else:
                self._embeddings = np.stack([o.full_embedding for o in outputs])

            # Build metadata from outputs
            self._epitope_metadata = {}
            for output in outputs:
                self._epitope_metadata[output.epitope_id] = {
                    'pdb_id': output.pdb_id,
                    'antigen_chains': list(output.chain_embeddings.keys()),
                    'epitope_residues': {
                        chain_id: result.epitope_indices
                        for chain_id, result in output.chain_embeddings.items()
                    },
                    'total_residues': output.total_epitope_residues
                }

            logger.info(f"Built index with {len(outputs)} encoder outputs")

        except Exception as e:
            raise GroupingError(f"Failed to build index: {e}") from e

    def compute_similarity_matrix(self) -> np.ndarray:
        """
        Compute pairwise cosine similarity matrix.

        For L2-normalized vectors: cosine_sim = dot_product

        Returns:
            Similarity matrix (N, N)
        """
        if self._embeddings is None:
            raise GroupingError("Index not built. Call build_index() first.")

        # For L2-normalized vectors, cosine similarity = dot product
        self._similarity_matrix = self._embeddings @ self._embeddings.T

        return self._similarity_matrix

    def find_groups(
        self,
        similarity_threshold: Optional[float] = None,
        antigen_clusters: Optional[Dict[str, str]] = None
    ) -> List[EpitopeGroup]:
        """
        Find connected components of similar epitopes.

        Args:
            similarity_threshold: Override default threshold
            antigen_clusters: Optional mapping epitope_id -> cluster_id for filtering

        Returns:
            List of EpitopeGroup objects
        """
        if self._embeddings is None:
            raise GroupingError("Index not built. Call build_index() first.")

        threshold = similarity_threshold or self.similarity_threshold

        # Compute similarity matrix if not already done
        if self._similarity_matrix is None:
            self.compute_similarity_matrix()

        # Build graph
        G = nx.Graph()
        n = len(self._epitope_ids)

        # Add all nodes
        for i, epitope_id in enumerate(self._epitope_ids):
            G.add_node(epitope_id, index=i, pdb_id=self._pdb_ids[i])

        # Add edges based on similarity
        for i in range(n):
            for j in range(i + 1, n):
                sim = self._similarity_matrix[i, j]

                # Check threshold
                if sim < threshold:
                    continue

                # Exclude same PDB
                if self.exclude_same_pdb and self._pdb_ids[i] == self._pdb_ids[j]:
                    continue

                # Exclude same antigen cluster if provided
                if antigen_clusters:
                    cluster_i = antigen_clusters.get(self._epitope_ids[i])
                    cluster_j = antigen_clusters.get(self._epitope_ids[j])
                    if cluster_i and cluster_j and cluster_i == cluster_j:
                        continue

                G.add_edge(self._epitope_ids[i], self._epitope_ids[j], similarity=sim)

        # Find connected components
        components = list(nx.connected_components(G))

        # Filter by minimum size and create EpitopeGroup objects
        groups = []
        group_idx = 0

        for component in components:
            if len(component) < self.min_group_size:
                continue

            members = list(component)

            # Select reference: node with highest degree (most connections)
            subgraph = G.subgraph(members)
            degrees = dict(subgraph.degree())
            reference_id = max(degrees, key=degrees.get)

            # Calculate average similarity within group
            similarities = []
            for u, v, data in subgraph.edges(data=True):
                similarities.append(data['similarity'])
            avg_sim = np.mean(similarities) if similarities else 0.0

            group = EpitopeGroup(
                group_id=f"group_{group_idx:04d}",
                reference_epitope_id=reference_id,
                member_epitope_ids=members,
                avg_similarity=float(avg_sim),
                metadata={'num_edges': len(subgraph.edges())}
            )
            groups.append(group)
            group_idx += 1

        logger.info(f"Found {len(groups)} groups from {n} epitopes")

        return groups

    def find_groups_detailed(
        self,
        similarity_threshold: Optional[float] = None
    ) -> GroupingOutput:
        """
        Find groups with detailed output including member info.

        Returns:
            GroupingOutput with full group details
        """
        if self._embeddings is None:
            raise GroupingError("Index not built. Call build_index() first.")

        threshold = similarity_threshold or self.similarity_threshold

        # Compute similarity matrix if not already done
        if self._similarity_matrix is None:
            self.compute_similarity_matrix()

        # Build graph
        G = nx.Graph()
        n = len(self._epitope_ids)

        # Create index lookup
        id_to_idx = {eid: i for i, eid in enumerate(self._epitope_ids)}

        # Add all nodes
        for i, epitope_id in enumerate(self._epitope_ids):
            G.add_node(epitope_id, index=i, pdb_id=self._pdb_ids[i])

        # Add edges based on similarity
        for i in range(n):
            for j in range(i + 1, n):
                sim = self._similarity_matrix[i, j]

                if sim < threshold:
                    continue

                if self.exclude_same_pdb and self._pdb_ids[i] == self._pdb_ids[j]:
                    continue

                G.add_edge(self._epitope_ids[i], self._epitope_ids[j], similarity=sim)

        # Find connected components
        components = list(nx.connected_components(G))

        # Build detailed groups
        groups = []
        group_idx = 0
        grouped_epitopes = set()

        for component in components:
            if len(component) < self.min_group_size:
                continue

            members_list = list(component)
            grouped_epitopes.update(members_list)

            # Select reference
            subgraph = G.subgraph(members_list)
            degrees = dict(subgraph.degree())
            reference_id = max(degrees, key=degrees.get)
            ref_idx = id_to_idx[reference_id]

            # Calculate similarities
            similarities = []
            for u, v, data in subgraph.edges(data=True):
                similarities.append(data['similarity'])

            avg_sim = np.mean(similarities) if similarities else 0.0
            min_sim = np.min(similarities) if similarities else 0.0
            max_sim = np.max(similarities) if similarities else 0.0

            # Build member list
            members = []
            for epitope_id in members_list:
                idx = id_to_idx[epitope_id]
                is_ref = (epitope_id == reference_id)

                # Get similarity to reference
                sim_to_ref = None
                if not is_ref:
                    sim_to_ref = float(self._similarity_matrix[idx, ref_idx])

                # Get metadata
                meta = self._epitope_metadata.get(epitope_id, {})

                member = GroupMember(
                    epitope_id=epitope_id,
                    pdb_id=self._pdb_ids[idx],
                    is_reference=is_ref,
                    antigen_chains=meta.get('antigen_chains', []),
                    epitope_residues=meta.get('epitope_residues', {}),
                    total_residues=meta.get('total_residues', 0),
                    similarity_to_ref=sim_to_ref
                )
                members.append(member)

            # Sort: reference first, then by similarity
            members.sort(key=lambda m: (not m.is_reference, -(m.similarity_to_ref or 1.0)))

            group = GroupResult(
                group_id=f"group_{group_idx:04d}",
                reference_epitope_id=reference_id,
                member_count=len(members),
                avg_similarity=float(avg_sim),
                min_similarity=float(min_sim),
                max_similarity=float(max_sim),
                members=members
            )
            groups.append(group)
            group_idx += 1

        # Count singletons
        singleton_count = n - len(grouped_epitopes)

        return GroupingOutput(
            groups=groups,
            total_epitopes=n,
            grouped_epitopes=len(grouped_epitopes),
            singleton_count=singleton_count,
            similarity_threshold=threshold
        )

    def save_similarity_matrix(self, output_path: Path) -> None:
        """
        Save similarity matrix to HDF5.

        Args:
            output_path: Output HDF5 file path
        """
        if self._similarity_matrix is None:
            self.compute_similarity_matrix()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(output_path, 'w') as f:
            # Save full matrix
            f.create_dataset(
                'similarity_matrix',
                data=self._similarity_matrix,
                compression='gzip',
                compression_opts=4
            )

            # Save epitope IDs
            f.create_dataset(
                'epitope_ids',
                data=np.array(self._epitope_ids, dtype='S')
            )

            # Save PDB IDs
            f.create_dataset(
                'pdb_ids',
                data=np.array(self._pdb_ids, dtype='S')
            )

            # Metadata
            f.attrs['num_epitopes'] = len(self._epitope_ids)
            f.attrs['embedding_dim'] = self._embeddings.shape[1]

        logger.info(f"Saved similarity matrix to {output_path}")

    def save_sparse_similarity(
        self,
        output_path: Path,
        threshold: Optional[float] = None
    ) -> None:
        """
        Save sparse similarity matrix (only values above threshold).

        Args:
            output_path: Output HDF5 file path
            threshold: Minimum similarity to save (default: self.similarity_threshold)
        """
        if self._similarity_matrix is None:
            self.compute_similarity_matrix()

        threshold = threshold or self.similarity_threshold
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Find pairs above threshold
        n = len(self._epitope_ids)
        rows, cols, values = [], [], []

        for i in range(n):
            for j in range(i + 1, n):
                sim = self._similarity_matrix[i, j]
                if sim >= threshold:
                    rows.append(i)
                    cols.append(j)
                    values.append(sim)

        with h5py.File(output_path, 'w') as f:
            f.create_dataset('rows', data=np.array(rows, dtype=np.int32))
            f.create_dataset('cols', data=np.array(cols, dtype=np.int32))
            f.create_dataset('values', data=np.array(values, dtype=np.float32))
            f.create_dataset('epitope_ids', data=np.array(self._epitope_ids, dtype='S'))
            f.create_dataset('pdb_ids', data=np.array(self._pdb_ids, dtype='S'))

            f.attrs['num_epitopes'] = n
            f.attrs['num_pairs'] = len(values)
            f.attrs['threshold'] = threshold

        logger.info(f"Saved sparse similarity ({len(values)} pairs) to {output_path}")


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_groups_json(output: GroupingOutput, output_path: Path) -> None:
    """
    Save grouping output to JSON file.

    Args:
        output: GroupingOutput object
        output_path: Output JSON file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(output.to_dict(), f, indent=2, cls=NumpyEncoder)

    logger.info(f"Saved {len(output.groups)} groups to {output_path}")


def save_grouping_stats_csv(output: GroupingOutput, output_path: Path) -> None:
    """
    Save grouping statistics to CSV.

    Args:
        output: GroupingOutput object
        output_path: Output CSV file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'group_id', 'member_count', 'avg_similarity',
            'min_similarity', 'max_similarity', 'reference_pdb', 'reference_epitope'
        ])

        for group in output.groups:
            ref_pdb = next(
                (m.pdb_id for m in group.members if m.is_reference),
                ''
            )
            writer.writerow([
                group.group_id,
                group.member_count,
                f"{group.avg_similarity:.4f}",
                f"{group.min_similarity:.4f}",
                f"{group.max_similarity:.4f}",
                ref_pdb,
                group.reference_epitope_id
            ])

    logger.info(f"Saved grouping stats to {output_path}")


def generate_grouping_report(output: GroupingOutput) -> str:
    """
    Generate human-readable grouping report.

    Args:
        output: GroupingOutput object

    Returns:
        Formatted report string
    """
    lines = []
    lines.append("=" * 70)
    lines.append("Epitope Grouping Report")
    lines.append("=" * 70)

    lines.append(f"\nSimilarity threshold: {output.similarity_threshold}")
    lines.append(f"Total epitopes: {output.total_epitopes}")
    lines.append(f"Grouped epitopes: {output.grouped_epitopes} ({output.grouped_epitopes/output.total_epitopes*100:.1f}%)")
    lines.append(f"Singletons (no matches): {output.singleton_count} ({output.singleton_count/output.total_epitopes*100:.1f}%)")
    lines.append(f"Number of groups: {len(output.groups)}")

    if output.groups:
        group_sizes = [g.member_count for g in output.groups]
        avg_similarities = [g.avg_similarity for g in output.groups]

        lines.append(f"\nGroup size statistics:")
        lines.append(f"  Mean: {np.mean(group_sizes):.1f}")
        lines.append(f"  Median: {np.median(group_sizes):.1f}")
        lines.append(f"  Range: {min(group_sizes)} - {max(group_sizes)}")

        lines.append(f"\nSimilarity statistics:")
        lines.append(f"  Mean avg_similarity: {np.mean(avg_similarities):.4f}")
        lines.append(f"  Range: {min(avg_similarities):.4f} - {max(avg_similarities):.4f}")

        lines.append(f"\nTop 5 largest groups:")
        sorted_groups = sorted(output.groups, key=lambda g: -g.member_count)[:5]
        for g in sorted_groups:
            lines.append(f"  {g.group_id}: {g.member_count} members, avg_sim={g.avg_similarity:.4f}")

    lines.append("=" * 70)

    return '\n'.join(lines)
