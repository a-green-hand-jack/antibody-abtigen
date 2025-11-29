"""
UniProt and Ensembl API mapping module.

This module handles:
1. PDB chain to UniProt ID mapping via SIFTS
2. UniProt ID to ortholog mapping via Ensembl
3. Sequence identity calculation between human and mouse proteins
"""

import os
import time
import requests
from typing import Optional, Dict, List, Tuple
from functools import lru_cache


# API endpoints
SIFTS_MAPPING_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}"
ENSEMBL_HOMOLOGY_URL = "https://rest.ensembl.org/homology/id/{ensembl_id}"
UNIPROT_IDMAPPING_URL = "https://rest.uniprot.org/idmapping"

# Taxonomy IDs
HUMAN_TAXON = 9606
MOUSE_TAXON = 10090


class RateLimiter:
    """Simple rate limiter to avoid API throttling."""

    def __init__(self, requests_per_second: float = 3.0):
        self.min_interval = 1.0 / requests_per_second
        self.last_request_time = 0

    def wait(self):
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request_time = time.time()


# Global rate limiters for different APIs
uniprot_limiter = RateLimiter(3.0)
ensembl_limiter = RateLimiter(10.0)
sifts_limiter = RateLimiter(5.0)


@lru_cache(maxsize=1000)
def get_uniprot_from_pdb_chain(pdb_id: str, chain_id: str) -> Optional[str]:
    """
    Get UniProt accession from PDB ID and chain ID using SIFTS.

    Args:
        pdb_id: PDB ID (e.g., "5FUO")
        chain_id: Chain ID (e.g., "A")

    Returns:
        UniProt accession or None if not found
    """
    sifts_limiter.wait()

    url = SIFTS_MAPPING_URL.format(pdb_id=pdb_id.lower())

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        data = response.json()

        # Parse SIFTS response
        pdb_data = data.get(pdb_id.lower(), {})
        uniprot_mappings = pdb_data.get('UniProt', {})

        for uniprot_id, mapping_info in uniprot_mappings.items():
            mappings = mapping_info.get('mappings', [])
            for mapping in mappings:
                if mapping.get('chain_id') == chain_id:
                    return uniprot_id

        return None

    except requests.RequestException as e:
        print(f"Error fetching SIFTS mapping for {pdb_id}_{chain_id}: {e}")
        return None


@lru_cache(maxsize=1000)
def get_uniprot_info(accession: str) -> Optional[Dict]:
    """
    Get detailed UniProt entry information.

    Args:
        accession: UniProt accession

    Returns:
        Dictionary with UniProt info or None
    """
    uniprot_limiter.wait()

    url = f"{UNIPROT_ENTRY_URL.format(accession=accession)}.json"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 404:
            return None
        response.raise_for_status()
        data = response.json()

        # Extract key information
        genes = data.get('genes', [])
        gene_name = genes[0].get('geneName', {}).get('value') if genes else None

        organism = data.get('organism', {})
        taxon_id = organism.get('taxonId')

        sequence = data.get('sequence', {})
        seq_value = sequence.get('value', '')

        # Get cross-references for Ensembl
        xrefs = data.get('uniProtKBCrossReferences', [])
        ensembl_ids = []
        for xref in xrefs:
            if xref.get('database') == 'Ensembl':
                ensembl_ids.append(xref.get('id'))

        return {
            'accession': accession,
            'gene_name': gene_name,
            'taxon_id': taxon_id,
            'organism': organism.get('scientificName'),
            'sequence': seq_value,
            'sequence_length': len(seq_value),
            'ensembl_ids': ensembl_ids,
            'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
        }

    except requests.RequestException as e:
        print(f"Error fetching UniProt info for {accession}: {e}")
        return None


@lru_cache(maxsize=500)
def find_mouse_ortholog(human_uniprot_id: str) -> Optional[Dict]:
    """
    Find mouse ortholog for a human UniProt entry.

    Uses UniProt's ortholog information first, then falls back to Ensembl.

    Args:
        human_uniprot_id: Human UniProt accession

    Returns:
        Dictionary with mouse ortholog info or None
    """
    # First, get the human protein info
    human_info = get_uniprot_info(human_uniprot_id)
    if not human_info or human_info.get('taxon_id') != HUMAN_TAXON:
        return None

    gene_name = human_info.get('gene_name')
    if not gene_name:
        return None

    # Search for mouse ortholog by gene name in UniProt
    uniprot_limiter.wait()

    query = f'(gene_exact:{gene_name}) AND (organism_id:{MOUSE_TAXON}) AND (reviewed:true)'
    params = {
        'query': query,
        'format': 'json',
        'size': 5,
        'fields': 'accession,gene_names,organism_name,sequence'
    }

    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        results = data.get('results', [])
        if not results:
            # Try alternative: search by gene name without exact match
            return _search_mouse_ortholog_fallback(gene_name, human_info)

        # Get the first reviewed mouse entry
        for result in results:
            mouse_acc = result.get('primaryAccession')
            mouse_info = get_uniprot_info(mouse_acc)
            if mouse_info:
                # Calculate sequence identity
                identity = calculate_sequence_identity(
                    human_info.get('sequence', ''),
                    mouse_info.get('sequence', '')
                )
                mouse_info['sequence_identity'] = identity
                mouse_info['human_ortholog'] = human_uniprot_id
                return mouse_info

        return None

    except requests.RequestException as e:
        print(f"Error searching for mouse ortholog of {human_uniprot_id}: {e}")
        return None


def _search_mouse_ortholog_fallback(gene_name: str, human_info: Dict) -> Optional[Dict]:
    """
    Fallback search for mouse ortholog using broader criteria.
    """
    uniprot_limiter.wait()

    # Try with gene name contains instead of exact
    query = f'(gene:{gene_name}) AND (organism_id:{MOUSE_TAXON}) AND (reviewed:true)'
    params = {
        'query': query,
        'format': 'json',
        'size': 10,
        'fields': 'accession,gene_names,organism_name,sequence'
    }

    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        results = data.get('results', [])

        best_match = None
        best_identity = 0

        for result in results:
            mouse_acc = result.get('primaryAccession')
            mouse_info = get_uniprot_info(mouse_acc)
            if mouse_info:
                identity = calculate_sequence_identity(
                    human_info.get('sequence', ''),
                    mouse_info.get('sequence', '')
                )
                if identity > best_identity:
                    best_identity = identity
                    best_match = mouse_info
                    best_match['sequence_identity'] = identity
                    best_match['human_ortholog'] = human_info['accession']

        return best_match

    except requests.RequestException:
        return None


def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """
    Calculate sequence identity between two sequences using simple alignment.

    For a more accurate calculation, we use a simple global alignment approach.
    For production, consider using Biopython's pairwise aligner.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Sequence identity as percentage (0-100)
    """
    if not seq1 or not seq2:
        return 0.0

    # Simple approach: count matching characters at aligned positions
    # For more accuracy, we should use proper alignment
    try:
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment

        # Use global alignment
        alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
        if not alignments:
            return 0.0

        alignment = alignments[0]
        aligned_seq1 = alignment[0]
        aligned_seq2 = alignment[1]

        # Count matches
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        alignment_length = len(aligned_seq1)

        if alignment_length == 0:
            return 0.0

        return (matches / alignment_length) * 100

    except ImportError:
        # Fallback: simple identity calculation
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        return (matches / max_len) * 100 if max_len > 0 else 0.0


@lru_cache(maxsize=500)
def find_pdb_structures_for_uniprot(uniprot_id: str) -> List[Dict]:
    """
    Find PDB structures containing a given UniProt protein.

    Args:
        uniprot_id: UniProt accession

    Returns:
        List of PDB entries with chain information
    """
    uniprot_limiter.wait()

    # Use RCSB PDB API to search for structures
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "polymer_entity",
        "request_options": {
            "return_all_hits": True
        }
    }

    try:
        response = requests.post(url, json=query, timeout=30)
        if response.status_code == 204:  # No content
            return []
        response.raise_for_status()
        data = response.json()

        results = []
        for result in data.get('result_set', []):
            entity_id = result.get('identifier', '')
            if '_' in entity_id:
                pdb_id, entity_num = entity_id.split('_')
                results.append({
                    'pdb_id': pdb_id,
                    'entity_id': entity_id,
                    'uniprot_id': uniprot_id
                })

        return results

    except requests.RequestException as e:
        print(f"Error searching PDB for {uniprot_id}: {e}")
        return []


def get_pdb_chain_for_entity(pdb_id: str, entity_id: str) -> List[str]:
    """
    Get chain IDs for a given PDB entity.

    Args:
        pdb_id: PDB ID
        entity_id: Entity ID (e.g., "1ABC_1")

    Returns:
        List of chain IDs
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id.split('_')[1]}"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 404:
            return []
        response.raise_for_status()
        data = response.json()

        # Get chain IDs from entity_poly
        entity_poly = data.get('entity_poly', {})
        chain_ids = entity_poly.get('pdbx_strand_id', '').split(',')
        return [c.strip() for c in chain_ids if c.strip()]

    except requests.RequestException:
        return []


if __name__ == "__main__":
    # Test the module
    print("Testing PDB to UniProt mapping...")
    uniprot_id = get_uniprot_from_pdb_chain("5FUO", "A")
    print(f"5FUO chain A -> UniProt: {uniprot_id}")

    if uniprot_id:
        print("\nTesting UniProt info retrieval...")
        info = get_uniprot_info(uniprot_id)
        if info:
            print(f"Gene: {info.get('gene_name')}")
            print(f"Organism: {info.get('organism')}")
            print(f"Protein: {info.get('protein_name')}")

        print("\nTesting mouse ortholog search...")
        mouse = find_mouse_ortholog(uniprot_id)
        if mouse:
            print(f"Mouse ortholog: {mouse.get('accession')}")
            print(f"Sequence identity: {mouse.get('sequence_identity'):.1f}%")
