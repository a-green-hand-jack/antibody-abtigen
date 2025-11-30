"""
Check chain IDs in original vs cleaned CIF files.
"""

from pathlib import Path
from Bio.PDB import MMCIFParser
import gemmi

def check_with_biopython(cif_path):
    """Check chains using Biopython."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("test", str(cif_path))
    model = structure[0]

    chains = {}
    for chain in model:
        chain_id = chain.get_id()
        num_residues = len([r for r in chain if r.get_id()[0] == ' '])
        chains[chain_id] = num_residues

    return chains

def check_with_gemmi(cif_path):
    """Check auth_asym_id and label_asym_id using Gemmi."""
    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()

    table = block.find_mmcif_category("_atom_site")
    if not table:
        return {}, {}

    tags = list(table.tags)

    try:
        idx_label = tags.index("_atom_site.label_asym_id")
        idx_auth = tags.index("_atom_site.auth_asym_id")
    except ValueError:
        return {}, {}

    label_chains = set()
    auth_chains = set()

    for row in table:
        label_chains.add(row[idx_label])
        auth_chains.add(row[idx_auth])

    return sorted(label_chains), sorted(auth_chains)

def main():
    print("=" * 70)
    print("Checking Chain IDs in CIF Files")
    print("=" * 70)

    # Check original
    original = Path("data/raw_cif/1a14.cif")
    print(f"\nOriginal: {original}")

    label_chains, auth_chains = check_with_gemmi(original)
    print(f"  label_asym_id: {label_chains}")
    print(f"  auth_asym_id: {auth_chains}")

    chains = check_with_biopython(original)
    print(f"  Biopython chains: {list(chains.keys())}")
    for chain_id, num_res in chains.items():
        print(f"    {chain_id}: {num_res} residues")

    # Check cleaned
    cleaned = Path("data/epitope_pipeline/cleaned/1A14_cleaned.cif")
    if cleaned.exists():
        print(f"\nCleaned: {cleaned}")

        label_chains, auth_chains = check_with_gemmi(cleaned)
        print(f"  label_asym_id: {label_chains}")
        print(f"  auth_asym_id: {auth_chains}")

        chains = check_with_biopython(cleaned)
        print(f"  Biopython chains: {list(chains.keys())}")
        for chain_id, num_res in chains.items():
            print(f"    {chain_id}: {num_res} residues")
    else:
        print(f"\nCleaned file not found: {cleaned}")

if __name__ == "__main__":
    main()
