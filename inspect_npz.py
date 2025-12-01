import numpy as np
import sys

path = "data/extract_epitope/pockets/6nf2_Q.npz"
try:
    data = np.load(path, allow_pickle=True)
    print(f"Keys: {list(data.keys())}")
    if 'res_seq_nums' in data:
        print(f"res_seq_nums: {data['res_seq_nums']}")
        print(f"Min: {np.min(data['res_seq_nums'])}, Max: {np.max(data['res_seq_nums'])}")
    if 'chain_ids' in data:
        print(f"chain_ids unique: {np.unique(data['chain_ids'])}")
except Exception as e:
    print(e)
