"""
Convert cleaned CIF files to PDB format for easier visualization.
"""

from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB format."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(cif_path))

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_path))

    print(f"✓ Converted: {cif_path.name} → {pdb_path.name}")

def main():
    cleaned_dir = Path("tests/cleaner/output_cleaned")
    pdb_dir = Path("tests/cleaner/output_pdb")
    pdb_dir.mkdir(parents=True, exist_ok=True)

    cif_files = sorted(cleaned_dir.glob("*_cleaned.cif"))

    print("=" * 70)
    print("Converting Cleaned CIF to PDB Format")
    print("=" * 70)
    print()

    for cif_path in cif_files:
        pdb_path = pdb_dir / cif_path.name.replace("_cleaned.cif", "_cleaned.pdb")
        try:
            convert_cif_to_pdb(cif_path, pdb_path)
        except Exception as e:
            print(f"✗ Failed: {cif_path.name} - {e}")

    print()
    print(f"PDB files saved to: {pdb_dir}")

if __name__ == "__main__":
    main()
