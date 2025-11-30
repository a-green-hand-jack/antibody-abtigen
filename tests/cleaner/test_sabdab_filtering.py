"""
Test how current cleaner handles different SAbDab antigen types.

This script checks:
1. Structures with antigen_chain=NA (no antigen)
2. Structures with antigen_type=peptide
3. Structures with antigen_type=hapten
"""

from pathlib import Path
import csv

def load_sabdab_summary():
    """Load SAbDab summary file."""
    sabdab_path = Path("data/meta/sabdab_summary_all.tsv")

    entries = {}
    with open(sabdab_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pdb_id = row['pdb'].upper()
            entries[pdb_id] = {
                'Hchain': row.get('Hchain', ''),
                'Lchain': row.get('Lchain', ''),
                'antigen_chain': row.get('antigen_chain', ''),
                'antigen_type': row.get('antigen_type', ''),
            }

    return entries

def main():
    print("=" * 70)
    print("SAbDab Antigen Type Analysis")
    print("=" * 70)

    entries = load_sabdab_summary()

    # Check which structures we have in raw_cif/
    raw_cif_dir = Path("data/raw_cif")
    available_pdbs = {f.stem.upper() for f in raw_cif_dir.glob("*.cif")}

    print(f"\nTotal SAbDab entries: {len(entries)}")
    print(f"Available CIF files: {len(available_pdbs)}")

    # Categorize
    no_antigen = []
    peptide_antigen = []
    hapten_antigen = []
    protein_antigen = []

    for pdb_id, info in entries.items():
        if pdb_id not in available_pdbs:
            continue

        antigen_chain = info['antigen_chain']
        antigen_type = info['antigen_type']

        if antigen_chain == 'NA' or not antigen_chain:
            no_antigen.append(pdb_id)
        elif 'peptide' in antigen_type:
            peptide_antigen.append(pdb_id)
        elif 'hapten' in antigen_type:
            hapten_antigen.append(pdb_id)
        elif 'protein' in antigen_type:
            protein_antigen.append(pdb_id)

    print("\n" + "=" * 70)
    print("Available CIF Files by Antigen Type")
    print("=" * 70)

    print(f"\nNo antigen (antigen_chain=NA): {len(no_antigen)}")
    if no_antigen:
        print(f"  Examples: {', '.join(no_antigen[:5])}")

    print(f"\nPeptide antigen: {len(peptide_antigen)}")
    if peptide_antigen:
        print(f"  Examples: {', '.join(peptide_antigen[:5])}")
        # Show details for first example
        pdb_id = peptide_antigen[0]
        info = entries[pdb_id]
        print(f"  Details for {pdb_id}:")
        print(f"    Hchain: {info['Hchain']}")
        print(f"    Lchain: {info['Lchain']}")
        print(f"    antigen_chain: {info['antigen_chain']}")
        print(f"    antigen_type: {info['antigen_type']}")

    print(f"\nHapten antigen: {len(hapten_antigen)}")
    if hapten_antigen:
        print(f"  Examples: {', '.join(hapten_antigen[:5])}")

    print(f"\nProtein antigen: {len(protein_antigen)}")
    if protein_antigen:
        print(f"  Examples: {', '.join(protein_antigen[:5])}")

    print("\n" + "=" * 70)
    print("Current Behavior Analysis")
    print("=" * 70)

    print("\n1. No antigen (NA) cases:")
    print("   - SAbDab metadata: antigen_chain='' → split('|') = ['']")
    print("   - Current code: chain won't match '' → falls back to heuristics")
    print("   - Result: ⚠ May incorrectly classify some chains as 'antigen'")

    print("\n2. Peptide antigen cases:")
    print("   - SAbDab metadata: antigen_chain has value, antigen_type=peptide")
    print("   - Current code: Extracts peptide chain as 'antigen'")
    print("   - Result: ✓ Will include peptide, but may have very short sequence")
    print("   - Issue: Peptides may not be suitable for ESM-2 embedding")

    print("\n3. Hapten antigen cases:")
    print("   - SAbDab metadata: antigen_type=hapten (small molecules)")
    print("   - Current code: Will try to extract, but likely has no residues")
    print("   - Result: ⚠ May fail during sequence extraction")

if __name__ == "__main__":
    main()
