"""
Test loading cleaned CIF files to verify structure integrity.
"""

from pathlib import Path
from Bio.PDB import MMCIFParser
import numpy as np

def test_structure(cif_path):
    """Test if structure can be loaded and has valid coordinates."""
    print(f"\nTesting: {cif_path.name}")
    print(f"  File size: {cif_path.stat().st_size / 1024:.1f} KB")

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("test", str(cif_path))
        model = structure[0]

        total_atoms = 0
        total_residues = 0

        for chain in model:
            chain_id = chain.get_id()
            residues = [r for r in chain if r.get_id()[0] == ' ']
            atoms = list(chain.get_atoms())

            total_residues += len(residues)
            total_atoms += len(atoms)

            # Check coordinates
            coords = np.array([atom.get_coord() for atom in atoms])
            if len(coords) > 0:
                center = coords.mean(axis=0)
                print(f"  Chain {chain_id}: {len(residues)} residues, {len(atoms)} atoms")
                print(f"    Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")

        print(f"  Total: {total_residues} residues, {total_atoms} atoms")
        print(f"  ✓ Structure loaded successfully")
        return True

    except Exception as e:
        print(f"  ✗ Failed to load: {e}")
        return False

def main():
    print("=" * 70)
    print("Testing Cleaned CIF Files")
    print("=" * 70)

    cleaned_dir = Path("tests/cleaner/output_cleaned")
    cif_files = sorted(cleaned_dir.glob("*_cleaned.cif"))[:5]

    success_count = 0
    for cif_path in cif_files:
        if test_structure(cif_path):
            success_count += 1

    print("\n" + "=" * 70)
    print(f"Summary: {success_count}/{len(cif_files)} files loaded successfully")
    print("=" * 70)

if __name__ == "__main__":
    main()
