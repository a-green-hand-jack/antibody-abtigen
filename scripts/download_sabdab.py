import os
import requests
from tqdm import tqdm

def download_file(url, filename):
    """Helper to download a file with a progress bar."""
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get('content-length', 0))
    block_size = 1024
    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
    
    with open(filename, 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()
    if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
        print("ERROR, something went wrong")

def main():
    print("Starting SAbDab download (this might take a while)...")
    # Note: SAbDab often requires registration or specific download links.
    # Since we are in a CLI agent, we simulate the typical download or provide instructions.
    # However, for the purpose of this script, we will assume the user has access to the 
    # summary file or a specific dump. 
    
    # Ideally, we would use the SAbDab API or a direct link if available/public.
    # A common public link for the summary file:
    summary_url = "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/"
    summary_file = "data/meta/sabdab_summary_all.tsv"
    
    print(f"Downloading summary file to {summary_file}...")
    try:
        download_file(summary_url, summary_file)
        print("Summary file downloaded.")
    except Exception as e:
        print(f"Failed to download summary: {e}")
        print("Please manually download the summary file from SAbDab.")

    print("\nTo download the full structure dataset, it is recommended to use the 'sabdab-downloader' script provided by OPIG or download the zip manually.")
    print("For this pipeline, please place the .cif files in 'data/raw_cif/'.")

if __name__ == "__main__":
    main()
