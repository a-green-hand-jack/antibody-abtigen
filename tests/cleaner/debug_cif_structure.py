"""
Debug CIF structure issues by comparing original and cleaned files.
"""

from pathlib import Path
import gemmi

def analyze_cif_structure(cif_path):
    """Analyze CIF file structure details."""
    print(f"\nAnalyzing: {cif_path.name}")

    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()

    # Get atom_site table
    table = block.find_mmcif_category("_atom_site")
    if not table:
        print("  ✗ No _atom_site table found!")
        return

    tags = list(table.tags)
    print(f"  Total tags: {len(tags)}")
    print(f"  Total rows: {len(list(table))}")

    # Check important tags
    important_tags = [
        "_atom_site.group_PDB",
        "_atom_site.label_asym_id",
        "_atom_site.auth_asym_id",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]

    print("  Tag presence:")
    for tag in important_tags:
        present = tag in tags
        print(f"    {tag}: {'✓' if present else '✗'}")

    # Sample first few rows
    print("  First 3 atoms:")
    for i, row in enumerate(table):
        if i >= 3:
            break

        try:
            idx_x = tags.index("_atom_site.Cartn_x")
            idx_y = tags.index("_atom_site.Cartn_y")
            idx_z = tags.index("_atom_site.Cartn_z")
            idx_atom = tags.index("_atom_site.label_atom_id")
            idx_comp = tags.index("_atom_site.label_comp_id")
            idx_chain = tags.index("_atom_site.auth_asym_id")

            atom = row[idx_atom]
            comp = row[idx_comp]
            chain = row[idx_chain]
            x = row[idx_x]
            y = row[idx_y]
            z = row[idx_z]

            print(f"    {i+1}. {comp}:{atom} chain {chain} @ ({x}, {y}, {z})")
        except (ValueError, IndexError) as e:
            print(f"    {i+1}. Error reading row: {e}")

def main():
    print("=" * 70)
    print("Debugging CIF Structure")
    print("=" * 70)

    # Compare original and cleaned
    original = Path("data/raw_cif/1a14.cif")
    cleaned = Path("tests/cleaner/output_cleaned/1A14_cleaned.cif")

    analyze_cif_structure(original)
    analyze_cif_structure(cleaned)

if __name__ == "__main__":
    main()
