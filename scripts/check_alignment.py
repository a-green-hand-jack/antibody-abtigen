#!/usr/bin/env python3
"""
Check alignment quality using PyMOL.

This script loads the reference and aligned structures, and calculates:
1. Global RMSD (all CA atoms)
2. Epitope-only RMSD (CA atoms in epitope regions)
"""

import subprocess
import sys
import tempfile
from pathlib import Path
import json
import csv

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
ALIGNED_DIR = PROJECT_ROOT / "data/epitope_pipeline/aligned_test/group_0000"
ORCHESTRATOR_DIR = PROJECT_ROOT / "data/epitope_pipeline/orchestrator_test"

def load_epitope_residues():
    """Load epitope residues from CSV."""
    epitope_csv = ORCHESTRATOR_DIR / "embeddings/epitope_residues.csv"
    if not epitope_csv.exists():
        # Try other locations
        epitope_csv = PROJECT_ROOT / "data/epitope_pipeline/embeddings/csv/epitope_residues.csv"

    epitopes = {}
    with open(epitope_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = row['pdb_id']
            if pdb_id not in epitopes:
                epitopes[pdb_id] = []
            epitopes[pdb_id].append({
                'chain': row['chain_id'],
                'resi': row['auth_seq_id'],
                'resn': row['residue_type']
            })
    return epitopes

def create_pymol_script(reference_cif, aligned_cifs, epitopes, output_file):
    """Create PyMOL script to check alignment."""

    ref_pdb = Path(reference_cif).stem.split('_')[0]
    ref_epitope_residues = epitopes.get(ref_pdb, [])

    # Build reference epitope selection string
    if ref_epitope_residues:
        ref_epi_parts = [f"(chain {r['chain']} and resi {r['resi']})" for r in ref_epitope_residues]
        ref_epi_sel_base = " or ".join(ref_epi_parts)
        ref_epi_sel = f"ref and ({ref_epi_sel_base}) and name CA"
    else:
        ref_epi_sel = "none"

    script = f'''
# PyMOL script to check alignment quality
import json

results = {{}}

# Load reference
cmd.load("{reference_cif}", "ref")

print("Reference: {ref_pdb}")
print("Reference epitope residues: {len(ref_epitope_residues)}")

'''

    for aligned_cif in aligned_cifs:
        aligned_pdb = Path(aligned_cif).stem.split('_')[0]
        obj_name = f"aln_{aligned_pdb}"

        aligned_epitope = epitopes.get(aligned_pdb, [])

        # Build aligned epitope selection string
        if aligned_epitope:
            aln_epi_parts = [f"(chain {r['chain']} and resi {r['resi']})" for r in aligned_epitope]
            aln_epi_sel_base = " or ".join(aln_epi_parts)
            aln_epi_sel = f"{obj_name} and ({aln_epi_sel_base}) and name CA"
        else:
            aln_epi_sel = "none"

        script += f'''
# Load aligned structure: {aligned_pdb}
cmd.load("{aligned_cif}", "{obj_name}")

# Calculate global RMSD (all CA atoms) - align without moving
try:
    global_rmsd = cmd.align("ref and name CA", "{obj_name} and name CA", cycles=0)[0]
except Exception as e:
    print(f"Global align failed for {aligned_pdb}: {{e}}")
    global_rmsd = -1.0

# Calculate super RMSD (structure-based, sequence-independent)
try:
    super_result = cmd.super("{obj_name} and name CA", "ref and name CA", cycles=0)
    super_rmsd = super_result[0]
    super_atoms = super_result[1]
except Exception as e:
    print(f"Super failed for {aligned_pdb}: {{e}}")
    super_rmsd = -1.0
    super_atoms = 0

# Calculate epitope-only RMSD
try:
    epi_result = cmd.super("{aln_epi_sel}", "{ref_epi_sel}", cycles=0)
    epitope_rmsd = epi_result[0] if epi_result else -1.0
    epitope_atoms = epi_result[1] if epi_result else 0
except Exception as e:
    print(f"Epitope super failed for {aligned_pdb}: {{e}}")
    epitope_rmsd = -1.0
    epitope_atoms = 0

print(f"{aligned_pdb}: global_align={{global_rmsd:.3f}}, super={{super_rmsd:.3f}} ({{super_atoms}} atoms), epitope={{epitope_rmsd:.3f}} ({{epitope_atoms}} atoms)")

results["{aligned_pdb}"] = {{
    "global_align_rmsd": global_rmsd,
    "super_rmsd": super_rmsd,
    "super_atoms": super_atoms,
    "epitope_rmsd": epitope_rmsd,
    "epitope_atoms": epitope_atoms,
    "ref_epitope_count": {len(ref_epitope_residues)},
    "aligned_epitope_count": {len(aligned_epitope)}
}}

# Reload to reset position
cmd.delete("{obj_name}")
cmd.load("{aligned_cif}", "{obj_name}")

'''

    script += f'''
# Save results
with open("{output_file}", 'w') as f:
    json.dump(results, f, indent=2)

print("\\nResults saved to {output_file}")

# Visualize
cmd.show("cartoon", "all")
cmd.color("gray80", "ref")

# Color aligned structures
colors = ["red", "orange", "yellow", "green", "cyan", "violet", "magenta", "salmon", "lime"]
for i, obj in enumerate(cmd.get_object_list()[1:]):
    cmd.color(colors[i % len(colors)], obj)

cmd.zoom("all")
cmd.bg_color("white")
cmd.save("{ALIGNED_DIR}/alignment_check.pse")
print("\\nPyMOL session saved to {ALIGNED_DIR}/alignment_check.pse")
'''

    return script


def main():
    print("=" * 70)
    print("Checking Alignment Quality with PyMOL")
    print("=" * 70)

    # Load epitope residues
    epitopes = load_epitope_residues()
    print(f"\nLoaded epitope data for {len(epitopes)} structures")

    # Get reference and aligned structures
    ref_dir = ALIGNED_DIR / "reference"
    aln_dir = ALIGNED_DIR / "aligned"

    reference_antigen = list(ref_dir.glob("*_antigen.cif"))[0]
    aligned_antigens = sorted(aln_dir.glob("*_antigen.cif"))

    print(f"\nReference: {reference_antigen.name}")
    print(f"Aligned structures: {len(aligned_antigens)}")

    # Create output file path
    output_file = ALIGNED_DIR / "alignment_rmsd_results.json"

    # Create PyMOL script
    pymol_script = create_pymol_script(
        str(reference_antigen),
        [str(p) for p in aligned_antigens],
        epitopes,
        str(output_file)
    )

    # Write script to temp file
    script_file = ALIGNED_DIR / "check_alignment.pml"
    with open(script_file, 'w') as f:
        f.write(pymol_script)
    print(f"\nPyMOL script written to: {script_file}")

    # Run PyMOL
    print("\nRunning PyMOL...")
    result = subprocess.run(
        ["mamba", "run", "-n", "pymol-env", "python", "-c",
         f"import pymol; pymol.finish_launching(['pymol', '-cq']); exec(open('{script_file}').read())"],
        capture_output=True, text=True, cwd=str(PROJECT_ROOT)
    )

    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr)

    # Load and display results
    if output_file.exists():
        with open(output_file, 'r') as f:
            results = json.load(f)

        print("\n" + "=" * 70)
        print("RMSD Results Summary")
        print("=" * 70)
        print(f"{'PDB':<8} {'Global Align':<14} {'Super RMSD':<14} {'Super Atoms':<12} {'Epi RMSD':<12} {'Epi Atoms'}")
        print("-" * 70)

        for pdb, data in sorted(results.items()):
            global_rmsd = data['global_align_rmsd']
            super_rmsd = data['super_rmsd']
            super_atoms = data['super_atoms']
            epi_rmsd = data['epitope_rmsd']
            epi_atoms = data['epitope_atoms']

            print(f"{pdb:<8} {global_rmsd:>12.3f}Å {super_rmsd:>12.3f}Å {super_atoms:>10} {epi_rmsd:>10.3f}Å {epi_atoms:>8}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
