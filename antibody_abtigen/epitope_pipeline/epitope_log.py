"""
Epitope logging utilities for the epitope-centric pipeline.

This module provides CSV logging for:
1. Epitope residue definitions (which chains, which residues)
2. Embedding metadata (sequence lengths, epitope sizes)
3. Processing statistics

Output CSV format:
- epitope_residues.csv: Detailed per-residue information
- epitope_summary.csv: Per-structure summary
"""

import csv
import logging
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass

from .core import (
    CleanedStructure,
    EpitopeResidues,
    ChainMapping,
)
from .encoder import EncoderOutput

logger = logging.getLogger(__name__)


@dataclass
class EpitopeLogEntry:
    """A single log entry for epitope residues."""
    pdb_id: str
    epitope_id: str
    chain_id: str
    auth_seq_id: int  # PDB residue number
    zero_based_index: int  # ESM-2 index
    residue_type: str  # Single letter amino acid


def save_epitope_residues_csv(
    epitopes: List[EpitopeResidues],
    structures: Dict[str, CleanedStructure],
    output_path: Path
) -> None:
    """
    Save detailed epitope residue information to CSV.

    Each row represents one epitope residue with its chain, PDB number,
    and 0-based index.

    Args:
        epitopes: List of EpitopeResidues objects
        structures: Dict mapping pdb_id -> CleanedStructure
        output_path: Output CSV file path

    CSV columns:
        pdb_id, epitope_id, chain_id, auth_seq_id, zero_based_index, residue_type
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    entries = []

    for epitope in epitopes:
        structure = structures.get(epitope.pdb_id)
        if structure is None:
            logger.warning(f"Structure not found for {epitope.pdb_id}")
            continue

        for chain_id, auth_seq_ids in epitope.antigen_chains.items():
            # Find chain mapping
            chain_mapping = None
            for m in structure.chain_mappings:
                if m.original_chain_id == chain_id:
                    chain_mapping = m
                    break

            if chain_mapping is None:
                logger.warning(f"Chain mapping not found for {chain_id} in {epitope.pdb_id}")
                continue

            # Build inverse mapping: auth_seq_id -> 0-based index
            auth_to_idx = {v: k for k, v in chain_mapping.auth_seq_id_map.items()}

            for auth_seq_id in auth_seq_ids:
                zero_based_idx = auth_to_idx.get(auth_seq_id, -1)

                # Get residue type from sequence
                residue_type = 'X'
                if zero_based_idx >= 0 and zero_based_idx < len(chain_mapping.sequence):
                    residue_type = chain_mapping.sequence[zero_based_idx]

                entries.append(EpitopeLogEntry(
                    pdb_id=epitope.pdb_id,
                    epitope_id=epitope.epitope_id,
                    chain_id=chain_id,
                    auth_seq_id=auth_seq_id,
                    zero_based_index=zero_based_idx,
                    residue_type=residue_type
                ))

    # Write CSV
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'pdb_id', 'epitope_id', 'chain_id',
            'auth_seq_id', 'zero_based_index', 'residue_type'
        ])

        for entry in entries:
            writer.writerow([
                entry.pdb_id,
                entry.epitope_id,
                entry.chain_id,
                entry.auth_seq_id,
                entry.zero_based_index,
                entry.residue_type
            ])

    logger.info(f"Saved {len(entries)} epitope residue entries to {output_path}")


def save_epitope_summary_csv(
    epitopes: List[EpitopeResidues],
    structures: Dict[str, CleanedStructure],
    output_path: Path,
    encoder_outputs: Optional[List[EncoderOutput]] = None
) -> None:
    """
    Save per-structure epitope summary to CSV.

    Args:
        epitopes: List of EpitopeResidues objects
        structures: Dict mapping pdb_id -> CleanedStructure
        output_path: Output CSV file path
        encoder_outputs: Optional list of EncoderOutput for embedding metadata

    CSV columns:
        pdb_id, epitope_id, num_antigen_chains, num_antibody_chains,
        total_epitope_residues, antigen_chains, epitope_residues_by_chain,
        distance_threshold, embedding_generated
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build encoder output lookup
    encoder_lookup = {}
    if encoder_outputs:
        for output in encoder_outputs:
            encoder_lookup[output.epitope_id] = output

    rows = []
    for epitope in epitopes:
        structure = structures.get(epitope.pdb_id)

        # Count chains
        num_antigen_chains = len(epitope.antigen_chains)
        num_antibody_chains = len(epitope.antibody_chains)

        # Total residues
        total_residues = epitope.total_residue_count()

        # Format antigen chains
        antigen_chains_str = '|'.join(sorted(epitope.antigen_chains.keys()))

        # Format residues by chain
        residues_by_chain = []
        for chain_id in sorted(epitope.antigen_chains.keys()):
            residues = epitope.antigen_chains[chain_id]
            residues_by_chain.append(f"{chain_id}:{len(residues)}")
        residues_by_chain_str = '|'.join(residues_by_chain)

        # Check if embedding was generated
        embedding_generated = epitope.epitope_id in encoder_lookup

        # Get sequence lengths if structure available
        antigen_seq_lengths = []
        if structure:
            for m in structure.get_antigen_chains():
                antigen_seq_lengths.append(f"{m.original_chain_id}:{len(m.sequence)}")
        antigen_seq_lengths_str = '|'.join(antigen_seq_lengths)

        rows.append({
            'pdb_id': epitope.pdb_id,
            'epitope_id': epitope.epitope_id,
            'num_antigen_chains': num_antigen_chains,
            'num_antibody_chains': num_antibody_chains,
            'total_epitope_residues': total_residues,
            'num_contacts': epitope.num_contacts,
            'antigen_chains': antigen_chains_str,
            'epitope_residues_by_chain': residues_by_chain_str,
            'antigen_seq_lengths': antigen_seq_lengths_str,
            'distance_threshold': epitope.distance_threshold,
            'embedding_generated': 'Yes' if embedding_generated else 'No'
        })

    # Write CSV
    with open(output_path, 'w', newline='') as f:
        fieldnames = [
            'pdb_id', 'epitope_id', 'num_antigen_chains', 'num_antibody_chains',
            'total_epitope_residues', 'num_contacts', 'antigen_chains',
            'epitope_residues_by_chain', 'antigen_seq_lengths',
            'distance_threshold', 'embedding_generated'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Saved {len(rows)} epitope summaries to {output_path}")


def save_embedding_stats_csv(
    encoder_outputs: List[EncoderOutput],
    output_path: Path
) -> None:
    """
    Save embedding statistics to CSV.

    Args:
        encoder_outputs: List of EncoderOutput objects
        output_path: Output CSV file path

    CSV columns:
        pdb_id, epitope_id, num_chains, total_epitope_residues,
        chain_details, embedding_method
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for output in encoder_outputs:
        # Format chain details
        chain_details = []
        for chain_id, result in output.chain_embeddings.items():
            chain_details.append(
                f"{chain_id}:seq_len={len(result.sequence)},epitope={result.num_epitope_residues}"
            )
        chain_details_str = '|'.join(chain_details)

        rows.append({
            'pdb_id': output.pdb_id,
            'epitope_id': output.epitope_id,
            'num_chains': len(output.chain_embeddings),
            'total_epitope_residues': output.total_epitope_residues,
            'chain_details': chain_details_str,
            'embedding_method': output.embedding_method
        })

    # Write CSV
    with open(output_path, 'w', newline='') as f:
        fieldnames = [
            'pdb_id', 'epitope_id', 'num_chains', 'total_epitope_residues',
            'chain_details', 'embedding_method'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Saved {len(rows)} embedding stats to {output_path}")


def generate_epitope_report(
    epitopes: List[EpitopeResidues],
    encoder_outputs: Optional[List[EncoderOutput]] = None
) -> str:
    """
    Generate a human-readable summary report.

    Args:
        epitopes: List of EpitopeResidues
        encoder_outputs: Optional list of EncoderOutput

    Returns:
        Formatted report string
    """
    import statistics

    lines = []
    lines.append("=" * 70)
    lines.append("Epitope Extraction and Embedding Report")
    lines.append("=" * 70)

    # Epitope statistics
    total_epitopes = len(epitopes)
    residue_counts = [e.total_residue_count() for e in epitopes]
    contact_counts = [e.num_contacts for e in epitopes]

    lines.append(f"\nEpitope Statistics:")
    lines.append(f"  Total epitopes: {total_epitopes}")

    if residue_counts:
        lines.append(f"  Epitope residues:")
        lines.append(f"    Mean: {statistics.mean(residue_counts):.1f}")
        lines.append(f"    Median: {statistics.median(residue_counts):.1f}")
        lines.append(f"    Range: {min(residue_counts)} - {max(residue_counts)}")

    if contact_counts:
        lines.append(f"  Contacts:")
        lines.append(f"    Mean: {statistics.mean(contact_counts):.1f}")
        lines.append(f"    Range: {min(contact_counts)} - {max(contact_counts)}")

    # Multi-chain statistics
    multi_chain = sum(1 for e in epitopes if len(e.antigen_chains) > 1)
    lines.append(f"  Multi-chain epitopes: {multi_chain} ({multi_chain/total_epitopes*100:.1f}%)")

    # Embedding statistics
    if encoder_outputs:
        lines.append(f"\nEmbedding Statistics:")
        lines.append(f"  Total embeddings: {len(encoder_outputs)}")

        seq_lengths = []
        for output in encoder_outputs:
            for result in output.chain_embeddings.values():
                seq_lengths.append(len(result.sequence))

        if seq_lengths:
            lines.append(f"  Sequence lengths:")
            lines.append(f"    Mean: {statistics.mean(seq_lengths):.1f}")
            lines.append(f"    Range: {min(seq_lengths)} - {max(seq_lengths)}")

    lines.append("=" * 70)

    return '\n'.join(lines)
