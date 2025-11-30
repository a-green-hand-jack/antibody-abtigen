import pandas as pd
import requests
import os
from tqdm import tqdm
import concurrent.futures
import typer

def download_cif(pdb_id, output_dir):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    output_path = os.path.join(output_dir, f"{pdb_id}.cif")
    
    if os.path.exists(output_path):
        return "skipped"
        
    try:
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            return "downloaded"
        else:
            return f"failed_{response.status_code}"
    except Exception as e:
        return f"error_{str(e)}"

def main(
    summary_path: str = "data/meta/sabdab_summary_all.tsv",
    output_dir: str = "data/raw_cif",
    limit: int = -1,
    max_workers: int = 5
):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    print(f"Reading summary from {summary_path}...")
    df = pd.read_csv(summary_path, sep='\t')
    
    # Filter for entries with antigen chains (optional, but relevant for our task)
    # df = df[df['antigen_chain'] != 'NA']
    
    unique_pdbs = df['pdb'].unique()
    print(f"Found {len(unique_pdbs)} unique PDBs.")
    
    if limit > 0:
        unique_pdbs = unique_pdbs[:limit]
        print(f"Limiting to {limit} PDBs for testing.")
        
    print(f"Downloading to {output_dir} with {max_workers} workers...")
    
    results = {"skipped": 0, "downloaded": 0, "failed": 0}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Create a progress bar
        future_to_pdb = {executor.submit(download_cif, pdb, output_dir): pdb for pdb in unique_pdbs}
        
        for future in tqdm(concurrent.futures.as_completed(future_to_pdb), total=len(unique_pdbs)):
            pdb = future_to_pdb[future]
            try:
                status = future.result()
                if status == "skipped":
                    results["skipped"] += 1
                elif status == "downloaded":
                    results["downloaded"] += 1
                else:
                    results["failed"] += 1
                    # print(f"Failed {pdb}: {status}")
            except Exception as exc:
                print(f"{pdb} generated an exception: {exc}")
                results["failed"] += 1

    print("\nDownload Summary:")
    print(results)

if __name__ == "__main__":
    typer.run(main)
