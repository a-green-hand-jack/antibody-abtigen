import os
import argparse
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tqdm import tqdm
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def download_cif(pdb_id: str, output_dir: Path) -> bool:
    """
    Downloads the CIF file for a given PDB ID.

    Args:
        pdb_id: The PDB ID (4 characters).
        output_dir: The directory to save the file.

    Returns:
        True if downloaded or already exists, False if failed.
    """
    pdb_id = pdb_id.lower()
    filename = f"{pdb_id}.cif"
    file_path = output_dir / filename

    if file_path.exists():
        # logger.info(f"File {filename} already exists. Skipping.")
        return True

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(file_path, 'wb') as f:
                f.write(response.content)
            return True
        else:
            logger.error(f"Failed to download {pdb_id}. Status code: {response.status_code}")
            return False
    except Exception as e:
        logger.error(f"Error downloading {pdb_id}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Download CIF files for PDB IDs listed in a CSV file.")
    parser.add_argument(
        "input_csv",
        type=str,
        help="Path to the input CSV file containing a 'pdb' column."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Directory to save downloaded CIF files. Defaults to ./data/raw_data/<csv_filename_stem>/"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Number of parallel workers for downloading. Default is 16."
    )

    args = parser.parse_args()

    input_path = Path(args.input_csv)

    if not input_path.exists():
        logger.error(f"Input CSV file not found: {input_path}")
        return

    # Determine output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        # Default: ./data/raw_data/<stem>/
        # Note: Relative to current working directory, not script location, as is standard convention.
        output_dir = Path("data/raw_data") / input_path.stem

    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Read CSV
    try:
        df = pd.read_csv(input_path)
    except Exception as e:
        logger.error(f"Failed to read CSV file: {e}")
        return

    if 'pdb' not in df.columns:
        logger.error("The CSV file must contain a 'pdb' column.")
        return

    pdb_ids = df['pdb'].unique()
    logger.info(f"Found {len(pdb_ids)} unique PDB IDs to download.")

    # Download in parallel
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        # Create a list of tasks
        futures = [
            executor.submit(download_cif, pdb_id, output_dir)
            for pdb_id in pdb_ids
        ]

        # Monitor progress
        results = []
        for future in tqdm(futures, total=len(pdb_ids), desc="Downloading CIFs"):
            results.append(future.result())

    success_count = sum(results)
    logger.info(f"Download complete. Successfully processed {success_count}/{len(pdb_ids)} PDBs.")

if __name__ == "__main__":
    main()
