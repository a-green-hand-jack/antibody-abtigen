"""
ESM-2 epitope encoder module for the epitope-centric pipeline.

This module implements Phase 3: Epitope Embedding.

Key features:
1. Full-context embedding: Process entire antigen sequence through ESM-2
2. Multi-chain support: Encode each chain separately, weighted average by epitope residue count
3. Dual output: Both full-chain embeddings and epitope-specific embeddings
4. FP16 mixed precision for memory efficiency on A100

Embedding strategy:
1. For each antigen chain:
   - Pass full sequence through ESM-2 3B
   - Extract embeddings at epitope residue positions
   - Mean pool epitope embeddings
2. Weighted average across chains (by epitope residue count)
3. L2 normalize final embedding
"""

import logging
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field

import numpy as np
import torch
import esm

from .core import (
    EpitopeEncoder,
    EpitopeEmbedding,
    EpitopeResidues,
    CleanedStructure,
    ChainMapping,
    EmbeddingError,
)

logger = logging.getLogger(__name__)


@dataclass
class ChainEmbeddingResult:
    """Result of encoding a single antigen chain."""
    chain_id: str
    sequence: str
    full_embedding: np.ndarray  # (seq_len, 2560) or (2560,) if pooled
    epitope_indices: List[int]  # 0-based indices of epitope residues
    epitope_embedding: np.ndarray  # (2560,) mean-pooled epitope embedding
    num_epitope_residues: int


@dataclass
class EncoderOutput:
    """
    Complete output from the encoder.

    Contains both full-chain and epitope-specific embeddings for all chains.
    """
    epitope_id: str
    pdb_id: str

    # Per-chain results
    chain_embeddings: Dict[str, ChainEmbeddingResult]

    # Aggregated embeddings (weighted by epitope residue count)
    full_embedding: np.ndarray  # (2560,) weighted average of full chain embeddings
    epitope_embedding: np.ndarray  # (2560,) weighted average of epitope embeddings

    # Metadata
    total_epitope_residues: int
    embedding_method: str = "esm2_t36_3B_UR50D"


class ESM2EpitopeEncoder(EpitopeEncoder):
    """
    ESM-2 3B model encoder for epitope embeddings.

    Encodes each antigen chain separately and combines results weighted
    by the number of epitope residues in each chain.

    Usage:
        >>> encoder = ESM2EpitopeEncoder(device="cuda")
        >>> output = encoder.encode_full(epitope, structure)
        >>> print(output.epitope_embedding.shape)
        (2560,)
        >>> print(output.chain_embeddings.keys())
        dict_keys(['A', 'B'])
    """

    # ESM-2 3B model produces 2560-dimensional embeddings
    EMBEDDING_DIM = 2560

    def __init__(
        self,
        model_name: str = "esm2_t36_3B_UR50D",
        device: str = "cuda",
        use_fp16: bool = True,
        batch_size: int = 4,
        cache_dir: Optional[Path] = None
    ):
        """
        Initialize ESM-2 encoder.

        Args:
            model_name: ESM-2 model name (default: 3B parameter model)
            device: Device to run model on ('cuda' or 'cpu')
            use_fp16: Use FP16 mixed precision (recommended for A100)
            batch_size: Batch size for encoding multiple sequences
            cache_dir: Directory to cache model weights
        """
        self.model_name = model_name
        self.device = device
        self.use_fp16 = use_fp16 and device == "cuda"
        self.batch_size = batch_size
        self.cache_dir = cache_dir

        self._model = None
        self._alphabet = None
        self._batch_converter = None

    def _load_model(self):
        """Lazy load ESM-2 model."""
        if self._model is not None:
            return

        logger.info(f"Loading ESM-2 model: {self.model_name}")

        # Load model
        if self.cache_dir:
            torch.hub.set_dir(str(self.cache_dir))

        self._model, self._alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        self._batch_converter = self._alphabet.get_batch_converter()

        # Move to device
        self._model = self._model.to(self.device)
        self._model.eval()

        # Enable FP16 if requested
        if self.use_fp16:
            self._model = self._model.half()

        logger.info(f"ESM-2 model loaded on {self.device} (FP16: {self.use_fp16})")

    def _get_epitope_indices(
        self,
        chain_mapping: ChainMapping,
        epitope_auth_seq_ids: List[int]
    ) -> List[int]:
        """
        Convert PDB auth_seq_id to 0-based ESM-2 indices.

        Args:
            chain_mapping: Chain mapping with auth_seq_id_map
            epitope_auth_seq_ids: PDB residue numbers of epitope

        Returns:
            List of 0-based indices into the sequence
        """
        # Invert the mapping: auth_seq_id -> 0-based index
        auth_to_idx = {v: k for k, v in chain_mapping.auth_seq_id_map.items()}

        indices = []
        for auth_id in epitope_auth_seq_ids:
            if auth_id in auth_to_idx:
                indices.append(auth_to_idx[auth_id])
            else:
                logger.warning(
                    f"Epitope residue {auth_id} not found in chain {chain_mapping.original_chain_id} mapping"
                )

        return sorted(indices)

    def _encode_sequence(self, sequence: str) -> np.ndarray:
        """
        Encode a single sequence through ESM-2.

        Args:
            sequence: Amino acid sequence

        Returns:
            Per-residue embeddings (seq_len, 2560)
        """
        self._load_model()

        # Prepare batch (single sequence)
        data = [("seq", sequence)]
        batch_labels, batch_strs, batch_tokens = self._batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)

        # Forward pass
        with torch.no_grad():
            if self.use_fp16:
                with torch.amp.autocast('cuda'):
                    results = self._model(batch_tokens, repr_layers=[36])
            else:
                results = self._model(batch_tokens, repr_layers=[36])

        # Extract embeddings (remove BOS and EOS tokens)
        # Shape: (1, seq_len+2, 2560) -> (seq_len, 2560)
        embeddings = results["representations"][36][0, 1:-1, :].cpu().float().numpy()

        return embeddings

    def _encode_chain(
        self,
        chain_mapping: ChainMapping,
        epitope_auth_seq_ids: List[int]
    ) -> ChainEmbeddingResult:
        """
        Encode a single antigen chain.

        Args:
            chain_mapping: Chain mapping with sequence and index map
            epitope_auth_seq_ids: PDB residue numbers of epitope in this chain

        Returns:
            ChainEmbeddingResult with full and epitope embeddings
        """
        sequence = chain_mapping.sequence

        # Get 0-based indices for epitope residues
        epitope_indices = self._get_epitope_indices(chain_mapping, epitope_auth_seq_ids)

        if not epitope_indices:
            raise EmbeddingError(
                f"No valid epitope indices found for chain {chain_mapping.original_chain_id}"
            )

        # Encode full sequence
        full_embeddings = self._encode_sequence(sequence)  # (seq_len, 2560)

        # Extract epitope embeddings and mean pool
        epitope_vectors = full_embeddings[epitope_indices]  # (num_epitope, 2560)
        epitope_embedding = np.mean(epitope_vectors, axis=0)  # (2560,)

        # Mean pool full sequence for full embedding
        full_embedding = np.mean(full_embeddings, axis=0)  # (2560,)

        return ChainEmbeddingResult(
            chain_id=chain_mapping.original_chain_id,
            sequence=sequence,
            full_embedding=full_embedding,
            epitope_indices=epitope_indices,
            epitope_embedding=epitope_embedding,
            num_epitope_residues=len(epitope_indices)
        )

    def encode_full(
        self,
        epitope: EpitopeResidues,
        structure: CleanedStructure
    ) -> EncoderOutput:
        """
        Encode epitope with full output (per-chain and aggregated).

        This is the main encoding method that returns all embeddings.

        Args:
            epitope: Epitope residue definition
            structure: Cleaned structure with chain mappings

        Returns:
            EncoderOutput with all embeddings

        Raises:
            EmbeddingError: If encoding fails
        """
        try:
            chain_results = {}
            total_epitope_residues = 0

            # Encode each antigen chain that has epitope residues
            for chain_id, auth_seq_ids in epitope.antigen_chains.items():
                # Find chain mapping
                chain_mapping = None
                for m in structure.chain_mappings:
                    if m.original_chain_id == chain_id and m.chain_type == 'antigen':
                        chain_mapping = m
                        break

                if chain_mapping is None:
                    logger.warning(f"Chain mapping not found for antigen chain {chain_id}")
                    continue

                # Encode this chain
                result = self._encode_chain(chain_mapping, auth_seq_ids)
                chain_results[chain_id] = result
                total_epitope_residues += result.num_epitope_residues

            if not chain_results:
                raise EmbeddingError(f"No chains encoded for {epitope.pdb_id}")

            # Weighted average across chains (by epitope residue count)
            if len(chain_results) == 1:
                # Single chain - use directly
                single_result = list(chain_results.values())[0]
                aggregated_full = single_result.full_embedding
                aggregated_epitope = single_result.epitope_embedding
            else:
                # Multiple chains - weighted average
                weights = []
                full_embeddings = []
                epitope_embeddings = []

                for result in chain_results.values():
                    weights.append(result.num_epitope_residues)
                    full_embeddings.append(result.full_embedding)
                    epitope_embeddings.append(result.epitope_embedding)

                weights = np.array(weights, dtype=np.float32)
                weights = weights / weights.sum()

                full_embeddings = np.stack(full_embeddings)
                epitope_embeddings = np.stack(epitope_embeddings)

                aggregated_full = np.average(full_embeddings, axis=0, weights=weights)
                aggregated_epitope = np.average(epitope_embeddings, axis=0, weights=weights)

            # L2 normalize
            aggregated_full = aggregated_full / np.linalg.norm(aggregated_full)
            aggregated_epitope = aggregated_epitope / np.linalg.norm(aggregated_epitope)

            return EncoderOutput(
                epitope_id=epitope.epitope_id,
                pdb_id=epitope.pdb_id,
                chain_embeddings=chain_results,
                full_embedding=aggregated_full.astype(np.float32),
                epitope_embedding=aggregated_epitope.astype(np.float32),
                total_epitope_residues=total_epitope_residues,
                embedding_method=self.model_name
            )

        except Exception as e:
            raise EmbeddingError(f"Failed to encode {epitope.pdb_id}: {str(e)}") from e

    def encode(
        self,
        epitope: EpitopeResidues,
        structure: CleanedStructure
    ) -> EpitopeEmbedding:
        """
        Generate embedding for a single epitope (interface method).

        Returns only the epitope embedding as EpitopeEmbedding dataclass.
        Use encode_full() for complete output with per-chain embeddings.

        Args:
            epitope: Epitope residue definition
            structure: Cleaned structure

        Returns:
            EpitopeEmbedding with normalized vector
        """
        output = self.encode_full(epitope, structure)

        # Collect all epitope indices across chains
        all_indices = []
        for chain_result in output.chain_embeddings.values():
            all_indices.extend(chain_result.epitope_indices)

        return EpitopeEmbedding(
            epitope_id=output.epitope_id,
            pdb_id=output.pdb_id,
            embedding=output.epitope_embedding,
            full_sequence="",  # Not applicable for multi-chain
            epitope_residue_indices=sorted(all_indices),
            embedding_method=output.embedding_method,
            metadata={
                "total_epitope_residues": output.total_epitope_residues,
                "num_chains": len(output.chain_embeddings)
            }
        )

    def batch_encode(
        self,
        epitopes: List[Tuple[EpitopeResidues, CleanedStructure]]
    ) -> List[EpitopeEmbedding]:
        """
        Batch encode for GPU efficiency (interface method).

        Args:
            epitopes: List of (epitope, structure) pairs

        Returns:
            List of EpitopeEmbedding objects
        """
        results = []
        for epitope, structure in epitopes:
            try:
                embedding = self.encode(epitope, structure)
                results.append(embedding)
            except EmbeddingError as e:
                logger.warning(f"Skipping {epitope.pdb_id}: {e}")
                continue
        return results

    def batch_encode_full(
        self,
        epitopes: List[Tuple[EpitopeResidues, CleanedStructure]]
    ) -> List[EncoderOutput]:
        """
        Batch encode with full output.

        Args:
            epitopes: List of (epitope, structure) pairs

        Returns:
            List of EncoderOutput objects
        """
        results = []
        for epitope, structure in epitopes:
            try:
                output = self.encode_full(epitope, structure)
                results.append(output)
            except EmbeddingError as e:
                logger.warning(f"Skipping {epitope.pdb_id}: {e}")
                continue
        return results
