"""
HDF5 embedding storage module for the epitope-centric pipeline.

This module implements Phase 3 storage: Save embeddings for later retrieval.

Storage structure in HDF5:
/embeddings/
    /{epitope_id}/
        full_embedding: (2560,) float32 - L2 normalized full antigen embedding
        epitope_embedding: (2560,) float32 - L2 normalized epitope embedding
        /chains/{chain_id}/
            full_embedding: (2560,) float32 - Per-chain full embedding
            epitope_embedding: (2560,) float32 - Per-chain epitope embedding
            sequence: string attribute
            epitope_indices: (N,) int32 - 0-based indices

/metadata/
    /{epitope_id}: JSON string with full metadata

Features:
- Gzip compression for efficient storage
- Incremental writes (append mode)
- Fast retrieval by epitope_id
- Complete metadata preservation
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional, Iterator
from dataclasses import asdict

import numpy as np
import h5py

from .core import (
    EmbeddingStore,
    EpitopeEmbedding,
)
from .encoder import EncoderOutput, ChainEmbeddingResult

logger = logging.getLogger(__name__)


class HDF5EmbeddingStore(EmbeddingStore):
    """
    HDF5-based storage for epitope embeddings.

    Supports both simple EpitopeEmbedding (interface requirement) and
    full EncoderOutput (with per-chain embeddings).

    Usage:
        >>> store = HDF5EmbeddingStore()
        >>> store.save_encoder_outputs(outputs, Path("embeddings.h5"))
        >>> loaded = store.load_encoder_output(Path("embeddings.h5"), "1a14_epi")
    """

    # Compression settings
    COMPRESSION = "gzip"
    COMPRESSION_OPTS = 4

    def __init__(self, chunk_size: int = 100):
        """
        Initialize HDF5 store.

        Args:
            chunk_size: Number of embeddings per chunk for efficient I/O
        """
        self.chunk_size = chunk_size

    def save(self, embeddings: List[EpitopeEmbedding], output_path: Path) -> None:
        """
        Save EpitopeEmbedding objects to HDF5 (interface method).

        Args:
            embeddings: List of embeddings to save
            output_path: HDF5 file path
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(output_path, 'a') as f:
            # Create groups if needed
            if 'embeddings' not in f:
                f.create_group('embeddings')
            if 'metadata' not in f:
                f.create_group('metadata')

            for emb in embeddings:
                # Skip if already exists
                if emb.epitope_id in f['embeddings']:
                    logger.debug(f"Skipping existing: {emb.epitope_id}")
                    continue

                # Create epitope group
                grp = f['embeddings'].create_group(emb.epitope_id)

                # Save embedding
                grp.create_dataset(
                    'epitope_embedding',
                    data=emb.embedding,
                    compression=self.COMPRESSION,
                    compression_opts=self.COMPRESSION_OPTS
                )

                # Save attributes
                grp.attrs['pdb_id'] = emb.pdb_id
                grp.attrs['embedding_method'] = emb.embedding_method
                grp.attrs['epitope_indices'] = json.dumps(emb.epitope_residue_indices)

                # Save full metadata as JSON
                metadata = {
                    'epitope_id': emb.epitope_id,
                    'pdb_id': emb.pdb_id,
                    'full_sequence': emb.full_sequence,
                    'epitope_residue_indices': emb.epitope_residue_indices,
                    'embedding_method': emb.embedding_method,
                    'metadata': emb.metadata
                }
                f['metadata'].create_dataset(
                    emb.epitope_id,
                    data=json.dumps(metadata)
                )

        logger.info(f"Saved {len(embeddings)} embeddings to {output_path}")

    def load(self, input_path: Path) -> List[EpitopeEmbedding]:
        """
        Load all EpitopeEmbedding objects from HDF5 (interface method).

        Args:
            input_path: HDF5 file path

        Returns:
            List of loaded embeddings
        """
        input_path = Path(input_path)
        if not input_path.exists():
            raise FileNotFoundError(f"HDF5 file not found: {input_path}")

        embeddings = []
        with h5py.File(input_path, 'r') as f:
            if 'embeddings' not in f:
                return embeddings

            for epitope_id in f['embeddings'].keys():
                emb = self.get_by_id(input_path, epitope_id)
                if emb:
                    embeddings.append(emb)

        return embeddings

    def get_by_id(self, input_path: Path, epitope_id: str) -> Optional[EpitopeEmbedding]:
        """
        Load a single embedding by ID (interface method).

        Args:
            input_path: HDF5 file path
            epitope_id: ID of epitope to retrieve

        Returns:
            EpitopeEmbedding if found, None otherwise
        """
        input_path = Path(input_path)
        if not input_path.exists():
            return None

        with h5py.File(input_path, 'r') as f:
            if 'embeddings' not in f or epitope_id not in f['embeddings']:
                return None

            grp = f['embeddings'][epitope_id]

            # Load embedding
            embedding = grp['epitope_embedding'][:]

            # Load metadata
            if 'metadata' in f and epitope_id in f['metadata']:
                metadata = json.loads(f['metadata'][epitope_id][()])
            else:
                # Fallback to attributes
                metadata = {
                    'epitope_id': epitope_id,
                    'pdb_id': grp.attrs.get('pdb_id', ''),
                    'full_sequence': '',
                    'epitope_residue_indices': json.loads(grp.attrs.get('epitope_indices', '[]')),
                    'embedding_method': grp.attrs.get('embedding_method', ''),
                    'metadata': {}
                }

            return EpitopeEmbedding(
                epitope_id=metadata['epitope_id'],
                pdb_id=metadata['pdb_id'],
                embedding=embedding,
                full_sequence=metadata['full_sequence'],
                epitope_residue_indices=metadata['epitope_residue_indices'],
                embedding_method=metadata['embedding_method'],
                metadata=metadata.get('metadata', {})
            )

    def save_encoder_outputs(
        self,
        outputs: List[EncoderOutput],
        output_path: Path
    ) -> None:
        """
        Save full EncoderOutput objects with per-chain embeddings.

        This saves:
        - Aggregated full_embedding and epitope_embedding
        - Per-chain embeddings with sequences and indices

        Args:
            outputs: List of EncoderOutput objects
            output_path: HDF5 file path
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(output_path, 'a') as f:
            # Create groups if needed
            if 'embeddings' not in f:
                f.create_group('embeddings')
            if 'metadata' not in f:
                f.create_group('metadata')

            for output in outputs:
                # Skip if already exists
                if output.epitope_id in f['embeddings']:
                    logger.debug(f"Skipping existing: {output.epitope_id}")
                    continue

                # Create epitope group
                grp = f['embeddings'].create_group(output.epitope_id)

                # Save aggregated embeddings
                grp.create_dataset(
                    'full_embedding',
                    data=output.full_embedding,
                    compression=self.COMPRESSION,
                    compression_opts=self.COMPRESSION_OPTS
                )
                grp.create_dataset(
                    'epitope_embedding',
                    data=output.epitope_embedding,
                    compression=self.COMPRESSION,
                    compression_opts=self.COMPRESSION_OPTS
                )

                # Save attributes
                grp.attrs['pdb_id'] = output.pdb_id
                grp.attrs['embedding_method'] = output.embedding_method
                grp.attrs['total_epitope_residues'] = output.total_epitope_residues

                # Save per-chain embeddings
                chains_grp = grp.create_group('chains')
                for chain_id, chain_result in output.chain_embeddings.items():
                    chain_grp = chains_grp.create_group(chain_id)

                    chain_grp.create_dataset(
                        'full_embedding',
                        data=chain_result.full_embedding,
                        compression=self.COMPRESSION,
                        compression_opts=self.COMPRESSION_OPTS
                    )
                    chain_grp.create_dataset(
                        'epitope_embedding',
                        data=chain_result.epitope_embedding,
                        compression=self.COMPRESSION,
                        compression_opts=self.COMPRESSION_OPTS
                    )
                    chain_grp.create_dataset(
                        'epitope_indices',
                        data=np.array(chain_result.epitope_indices, dtype=np.int32)
                    )

                    chain_grp.attrs['sequence'] = chain_result.sequence
                    chain_grp.attrs['num_epitope_residues'] = chain_result.num_epitope_residues

                # Save full metadata as JSON
                metadata = {
                    'epitope_id': output.epitope_id,
                    'pdb_id': output.pdb_id,
                    'total_epitope_residues': output.total_epitope_residues,
                    'embedding_method': output.embedding_method,
                    'chains': {
                        chain_id: {
                            'sequence': result.sequence,
                            'epitope_indices': result.epitope_indices,
                            'num_epitope_residues': result.num_epitope_residues
                        }
                        for chain_id, result in output.chain_embeddings.items()
                    }
                }
                f['metadata'].create_dataset(
                    output.epitope_id,
                    data=json.dumps(metadata)
                )

        logger.info(f"Saved {len(outputs)} encoder outputs to {output_path}")

    def load_encoder_output(
        self,
        input_path: Path,
        epitope_id: str
    ) -> Optional[EncoderOutput]:
        """
        Load a single EncoderOutput by ID.

        Args:
            input_path: HDF5 file path
            epitope_id: ID of epitope to retrieve

        Returns:
            EncoderOutput if found, None otherwise
        """
        input_path = Path(input_path)
        if not input_path.exists():
            return None

        with h5py.File(input_path, 'r') as f:
            if 'embeddings' not in f or epitope_id not in f['embeddings']:
                return None

            grp = f['embeddings'][epitope_id]

            # Load aggregated embeddings
            full_embedding = grp['full_embedding'][:]
            epitope_embedding = grp['epitope_embedding'][:]

            # Load per-chain embeddings
            chain_embeddings = {}
            if 'chains' in grp:
                for chain_id in grp['chains'].keys():
                    chain_grp = grp['chains'][chain_id]
                    chain_embeddings[chain_id] = ChainEmbeddingResult(
                        chain_id=chain_id,
                        sequence=chain_grp.attrs['sequence'],
                        full_embedding=chain_grp['full_embedding'][:],
                        epitope_indices=chain_grp['epitope_indices'][:].tolist(),
                        epitope_embedding=chain_grp['epitope_embedding'][:],
                        num_epitope_residues=chain_grp.attrs['num_epitope_residues']
                    )

            return EncoderOutput(
                epitope_id=epitope_id,
                pdb_id=grp.attrs['pdb_id'],
                chain_embeddings=chain_embeddings,
                full_embedding=full_embedding,
                epitope_embedding=epitope_embedding,
                total_epitope_residues=grp.attrs['total_epitope_residues'],
                embedding_method=grp.attrs['embedding_method']
            )

    def list_epitope_ids(self, input_path: Path) -> List[str]:
        """
        List all epitope IDs in the HDF5 file.

        Args:
            input_path: HDF5 file path

        Returns:
            List of epitope IDs
        """
        input_path = Path(input_path)
        if not input_path.exists():
            return []

        with h5py.File(input_path, 'r') as f:
            if 'embeddings' not in f:
                return []
            return list(f['embeddings'].keys())

    def get_summary(self, input_path: Path) -> Dict:
        """
        Get summary statistics of stored embeddings.

        Args:
            input_path: HDF5 file path

        Returns:
            Summary dictionary
        """
        input_path = Path(input_path)
        if not input_path.exists():
            return {'error': 'File not found'}

        with h5py.File(input_path, 'r') as f:
            if 'embeddings' not in f:
                return {'count': 0}

            epitope_ids = list(f['embeddings'].keys())
            total_residues = []
            chain_counts = []

            for eid in epitope_ids:
                grp = f['embeddings'][eid]
                if 'total_epitope_residues' in grp.attrs:
                    total_residues.append(grp.attrs['total_epitope_residues'])
                if 'chains' in grp:
                    chain_counts.append(len(grp['chains']))

            return {
                'count': len(epitope_ids),
                'total_epitope_residues': {
                    'mean': np.mean(total_residues) if total_residues else 0,
                    'min': int(np.min(total_residues)) if total_residues else 0,
                    'max': int(np.max(total_residues)) if total_residues else 0,
                },
                'chains_per_epitope': {
                    'mean': np.mean(chain_counts) if chain_counts else 0,
                    'min': int(np.min(chain_counts)) if chain_counts else 0,
                    'max': int(np.max(chain_counts)) if chain_counts else 0,
                }
            }
