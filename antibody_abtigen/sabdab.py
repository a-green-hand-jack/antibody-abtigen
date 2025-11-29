"""
SAbDab (Structural Antibody Database) data download and parsing module.

This module handles:
1. Downloading the SAbDab summary file
2. Parsing and filtering entries based on antigen species and type
3. Extracting relevant chain information
"""

import os
import requests
import pandas as pd
from typing import Optional
from tqdm import tqdm


SABDAB_SUMMARY_URL = "https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/"


def download_sabdab_summary(output_dir: str, force_download: bool = False) -> str:
    """
    Download the SAbDab summary file.

    Args:
        output_dir: Directory to save the summary file
        force_download: If True, re-download even if file exists

    Returns:
        Path to the downloaded summary file
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "sabdab_summary.tsv")

    if os.path.exists(output_path) and not force_download:
        print(f"SAbDab summary file already exists at {output_path}")
        return output_path

    print("Downloading SAbDab summary file...")

    response = requests.get(SABDAB_SUMMARY_URL, stream=True)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))

    with open(output_path, 'wb') as f:
        with tqdm(total=total_size, unit='iB', unit_scale=True, desc="Downloading") as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))

    print(f"Downloaded SAbDab summary to {output_path}")
    return output_path


def parse_sabdab_summary(summary_path: str) -> pd.DataFrame:
    """
    Parse the SAbDab summary TSV file.

    Args:
        summary_path: Path to the summary file

    Returns:
        DataFrame with parsed data
    """
    df = pd.read_csv(summary_path, sep='\t', low_memory=False)
    print(f"Loaded {len(df)} entries from SAbDab summary")
    return df


def filter_human_antigen_complexes(
    df: pd.DataFrame,
    resolution_threshold: float = 2.5,
    require_paired: bool = True
) -> pd.DataFrame:
    """
    Filter SAbDab entries for human antigen protein complexes.

    Args:
        df: SAbDab summary DataFrame
        resolution_threshold: Maximum resolution in Angstroms (default: 2.5)
        require_paired: If True, only include entries with both heavy and light chains

    Returns:
        Filtered DataFrame
    """
    # Filter for human antigens
    mask = df['antigen_species'].str.contains('homo sapiens', case=False, na=False)
    filtered = df[mask].copy()
    print(f"After filtering for human antigens: {len(filtered)} entries")

    # Filter for protein antigens
    mask = filtered['antigen_type'].str.contains('protein', case=False, na=False)
    filtered = filtered[mask]
    print(f"After filtering for protein antigens: {len(filtered)} entries")

    # Filter by resolution
    if 'resolution' in filtered.columns:
        # Convert resolution to numeric, handling non-numeric values
        filtered['resolution_num'] = pd.to_numeric(filtered['resolution'], errors='coerce')
        mask = filtered['resolution_num'] <= resolution_threshold
        filtered = filtered[mask]
        print(f"After filtering for resolution <= {resolution_threshold}A: {len(filtered)} entries")

    # Filter for paired antibodies (both H and L chains present)
    if require_paired:
        mask = (filtered['Hchain'].notna() & (filtered['Hchain'] != 'NA') &
                filtered['Lchain'].notna() & (filtered['Lchain'] != 'NA'))
        filtered = filtered[mask]
        print(f"After filtering for paired H/L chains: {len(filtered)} entries")

    # Filter for entries with antigen chain information
    mask = filtered['antigen_chain'].notna() & (filtered['antigen_chain'] != 'NA')
    filtered = filtered[mask]
    print(f"After filtering for antigen chain info: {len(filtered)} entries")

    return filtered


def extract_chain_info(row: pd.Series) -> dict:
    """
    Extract chain information from a SAbDab entry.

    Args:
        row: A row from the SAbDab DataFrame

    Returns:
        Dictionary with chain information
    """
    return {
        'pdb_id': row['pdb'],
        'heavy_chain': row['Hchain'],
        'light_chain': row['Lchain'],
        'antigen_chain': row['antigen_chain'],
        'antigen_species': row['antigen_species'],
        'antigen_type': row['antigen_type'],
        'antigen_name': row.get('antigen_name', 'Unknown'),
        'resolution': row.get('resolution', 'N/A'),
    }


def get_unique_complexes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get unique antibody-antigen complexes, removing duplicates.

    Args:
        df: Filtered SAbDab DataFrame

    Returns:
        DataFrame with unique complexes
    """
    # Group by PDB ID and antigen chain to get unique complexes
    # Keep the first entry for each unique complex
    unique = df.drop_duplicates(subset=['pdb', 'antigen_chain'], keep='first')
    print(f"Unique antibody-antigen complexes: {len(unique)}")
    return unique


def save_filtered_summary(df: pd.DataFrame, output_path: str) -> str:
    """
    Save filtered SAbDab summary to CSV.

    Args:
        df: Filtered DataFrame
        output_path: Path to save the filtered CSV

    Returns:
        Path to the saved file
    """
    df.to_csv(output_path, index=False)
    print(f"Saved filtered summary to {output_path}")
    return output_path


def find_mouse_antigen_in_sabdab(
    df: pd.DataFrame,
    antigen_name: str,
    human_uniprot_id: Optional[str] = None,
    resolution_threshold: float = 2.5
) -> Optional[dict]:
    """
    Search for mouse antigen structures in SAbDab data.

    This function prioritizes experimental structures from SAbDab
    over predicted structures from AlphaFold/RCSB.

    Args:
        df: Full SAbDab DataFrame (including all species)
        antigen_name: Human antigen name to search for
        human_uniprot_id: Optional UniProt ID to help match orthologs
        resolution_threshold: Maximum resolution in Angstroms

    Returns:
        Best matching mouse structure info or None
    """
    if df is None or df.empty:
        return None

    # Normalize antigen name for matching
    antigen_name_lower = antigen_name.lower().strip()

    # Filter for mouse antigens
    mouse_mask = df['antigen_species'].str.contains('mus musculus|mouse', case=False, na=False)
    mouse_entries = df[mouse_mask].copy()

    if mouse_entries.empty:
        return None

    # Filter for protein antigens
    protein_mask = mouse_entries['antigen_type'].str.contains('protein', case=False, na=False)
    mouse_entries = mouse_entries[protein_mask]

    if mouse_entries.empty:
        return None

    # Try to match by antigen name (case-insensitive)
    # Only return a match if we find a name match - don't return arbitrary mouse structures
    if 'antigen_name' not in mouse_entries.columns or not antigen_name_lower:
        return None

    # Exact match first
    name_mask = mouse_entries['antigen_name'].str.lower().str.strip() == antigen_name_lower
    matched = mouse_entries[name_mask]

    # If no exact match, try partial match (antigen name contains search term)
    if matched.empty:
        name_mask = mouse_entries['antigen_name'].str.lower().str.contains(
            antigen_name_lower, na=False, regex=False
        )
        matched = mouse_entries[name_mask]

    # If still no match, return None - don't return arbitrary mouse structures
    if matched.empty:
        return None

    mouse_entries = matched

    # Filter by resolution
    if 'resolution' in mouse_entries.columns:
        mouse_entries['resolution_num'] = pd.to_numeric(mouse_entries['resolution'], errors='coerce')
        res_mask = mouse_entries['resolution_num'] <= resolution_threshold
        mouse_entries = mouse_entries[res_mask]

    if mouse_entries.empty:
        return None

    # Select best structure (lowest resolution)
    return select_best_mouse_structure(mouse_entries)


def select_best_mouse_structure(candidates: pd.DataFrame) -> Optional[dict]:
    """
    Select the best mouse antigen structure from candidates.

    Selection criteria:
    1. Lowest resolution (highest quality)
    2. Presence of both H and L chains (complete antibody)
    3. Has antigen chain information

    Args:
        candidates: DataFrame of candidate mouse structures

    Returns:
        Dictionary with best structure info or None
    """
    if candidates.empty:
        return None

    # Sort by resolution (ascending = best first)
    if 'resolution_num' in candidates.columns:
        sorted_candidates = candidates.sort_values('resolution_num', ascending=True)
    else:
        sorted_candidates = candidates

    # Prefer structures with complete antibody (H + L)
    complete_mask = (
        sorted_candidates['Hchain'].notna() & (sorted_candidates['Hchain'] != 'NA') &
        sorted_candidates['Lchain'].notna() & (sorted_candidates['Lchain'] != 'NA')
    )
    complete_structures = sorted_candidates[complete_mask]

    if not complete_structures.empty:
        best = complete_structures.iloc[0]
    else:
        best = sorted_candidates.iloc[0]

    return {
        'pdb': best['pdb'],
        'antigen_chain': best['antigen_chain'],
        'antigen_name': best.get('antigen_name', 'Unknown'),
        'antigen_species': best['antigen_species'],
        'resolution': best.get('resolution', 'N/A'),
        'Hchain': best.get('Hchain', 'NA'),
        'Lchain': best.get('Lchain', 'NA'),
        'source': 'SAbDab'
    }


if __name__ == "__main__":
    # Test the module
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

    summary_path = download_sabdab_summary(data_dir)
    df = parse_sabdab_summary(summary_path)

    print("\nColumn names:")
    print(df.columns.tolist())

    print("\nFiltering for human antigen complexes...")
    filtered = filter_human_antigen_complexes(df)

    print("\nSample entries:")
    print(filtered.head())
