import torch
import esm
import os

def main():
    print("Downloading ESM-2 3B model (esm2_t36_3B_UR50D)...")
    print("This will be cached in ~/.cache/torch/hub/checkpoints/ by default.")
    
    # This triggers the download
    try:
        model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        print("Download complete. Model loaded successfully.")
    except Exception as e:
        print(f"Error downloading model: {e}")

if __name__ == "__main__":
    main()
