import os
import argparse
import requests
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

def download_file(url, filename):
    """Helper to download a file."""
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        if os.path.exists(filename):
            os.remove(filename)
        return False

def download_structure(pdb_id, output_dir):
    """Downloads a Chothia-numbered PDB file from SAbDab."""
    filename = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    if os.path.exists(filename):
        return True, pdb_id  # Already exists

    # SAbDab URL for Chothia-numbered PDBs
    url = f"https://opig.stats.ox.ac.uk/webapps/abdb/entries/{pdb_id}/structure/{pdb_id}.pdb"
    success = download_file(url, filename)
    return success, pdb_id
def main():
    parser = argparse.ArgumentParser(description="Download SAbDab summary and Chothia-numbered PDB structures.")
    parser.add_argument("--summary-path", default="data/meta/sabdab_summary_all.tsv", help="Path to save/load the summary file")
    parser.add_argument("--pdb-dir", default="data/raw_data/chothia", help="Directory to save Chothia-numbered PDB files")
    parser.add_argument("--workers", type=int, default=16, help="Number of parallel downloads")
    args = parser.parse_args()

    # 1. Download Summary
    if not os.path.exists(args.summary_path):
        print(f"Downloading SAbDab summary to {args.summary_path}...")
        summary_url = "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/"
        if not download_file(summary_url, args.summary_path):
            print("Failed to download summary file.")
            return
    else:
        print(f"Using existing summary file: {args.summary_path}")

    # 2. Filter Data
    print("Filtering SAbDab data...")
    try:
        df = pd.read_csv(args.summary_path, sep='\t')
    except Exception as e:
        print(f"Error reading summary file: {e}")
        return

    # Filter for protein antigens and valid antigen chains
    filtered_df = df[
        (df['antigen_type'].astype(str).str.contains('protein', case=False, na=False)) &
        (df['antigen_chain'].notna()) &
        (df['antigen_chain'] != '')
    ]

    pdb_ids = filtered_df['pdb'].unique()
    print(f"Found {len(pdb_ids)} unique PDB IDs with protein antigens.")

    # 3. Download PDBs
    os.makedirs(args.pdb_dir, exist_ok=True)
    print(f"Downloading Chothia-numbered PDB files to {args.pdb_dir} with {args.workers} workers...")

    failed_ids = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(download_structure, pdb_id, args.pdb_dir): pdb_id for pdb_id in pdb_ids}

        with tqdm(total=len(pdb_ids)) as pbar:
            for future in as_completed(futures):
                success, pdb_id = future.result()
                if not success:
                    failed_ids.append(pdb_id)
                pbar.update(1)

    print(f"Download complete. {len(pdb_ids) - len(failed_ids)} success, {len(failed_ids)} failed.")
    if failed_ids:
        print(f"Failed PDBs: {failed_ids[:10]} ...")
        # Optionally save failed IDs to a file

if __name__ == "__main__":
    main()
