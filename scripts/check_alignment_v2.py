#!/usr/bin/env python3
"""
Check alignment quality using PyMOL - V2 with correct epitope selections.
Run with: mamba run -n pymol-env python check_alignment_v2.py
"""

import json
import sys
import csv
from pathlib import Path

# Initialize PyMOL in headless mode
import pymol
pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd

# Paths
PROJECT_ROOT = Path("/ibex/user/wuj0c/Projects/Protein/GRADIENT/antibody-abtigen")
base_dir = PROJECT_ROOT / "data/epitope_pipeline/aligned_test/group_0000"
epitope_csv = PROJECT_ROOT / "data/epitope_pipeline/orchestrator_test/embeddings/epitope_residues.csv"
output_file = base_dir / "alignment_rmsd_results_v2.json"

def load_epitope_residues():
    """Load epitope residues from CSV."""
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
            })
    return epitopes

def build_selection(pdb_id, obj_name, epitopes):
    """Build PyMOL selection string for epitope residues."""
    if pdb_id not in epitopes:
        return "none"

    residues = epitopes[pdb_id]
    # Group by chain
    by_chain = {}
    for r in residues:
        chain = r['chain']
        if chain not in by_chain:
            by_chain[chain] = []
        by_chain[chain].append(r['resi'])

    # Build selection
    parts = []
    for chain, resis in by_chain.items():
        resi_str = "+".join(resis)
        parts.append(f"(chain {chain} and resi {resi_str})")

    if parts:
        sel = " or ".join(parts)
        return f"{obj_name} and ({sel}) and name CA"
    return "none"

# Load epitope data
epitopes = load_epitope_residues()
print(f"Loaded epitope data for {len(epitopes)} structures")

# Reference
ref_cif = str(base_dir / "reference/1A14_antigen.cif")
cmd.load(ref_cif, "ref")
ref_epi_sel = build_selection("1A14", "ref", epitopes)
ref_epi_ca = cmd.count_atoms(ref_epi_sel)
print(f"\nReference: 1A14")
print(f"Reference epitope CA atoms: {ref_epi_ca}")
print(f"Reference selection: {ref_epi_sel[:100]}...")

# Aligned structures
aligned_dir = base_dir / "aligned"
aligned_files = sorted(aligned_dir.glob("*_antigen.cif"))

results = {}

print("\n" + "="*90)
print(f"{'PDB':<8} {'EpiCA':<8} {'Global(Align)':<14} {'Global(Super)':<14} {'SuperAtoms':<12} {'Epitope':<12} {'EpiAtoms'}")
print("-"*90)

for cif_path in aligned_files:
    pdb_id = cif_path.stem.split('_')[0]
    obj_name = f"aln_{pdb_id}"

    # Load aligned structure
    cmd.load(str(cif_path), obj_name)

    aln_epi_sel = build_selection(pdb_id, obj_name, epitopes)
    aln_epi_ca = cmd.count_atoms(aln_epi_sel)

    # Calculate global RMSD using align (sequence-based)
    try:
        align_result = cmd.align(f"{obj_name} and name CA", "ref and name CA", cycles=0)
        global_align_rmsd = align_result[0]
    except Exception as e:
        global_align_rmsd = -1.0

    # Reset position by reloading
    cmd.delete(obj_name)
    cmd.load(str(cif_path), obj_name)

    # Calculate global RMSD using super (structure-based, sequence-independent)
    try:
        super_result = cmd.super(f"{obj_name} and name CA", "ref and name CA", cycles=0)
        global_super_rmsd = super_result[0]
        super_atoms = super_result[1]
    except Exception as e:
        global_super_rmsd = -1.0
        super_atoms = 0

    # Reset position by reloading
    cmd.delete(obj_name)
    cmd.load(str(cif_path), obj_name)

    # Calculate epitope-only RMSD (structure-based on epitope atoms)
    try:
        if aln_epi_ca > 0 and ref_epi_ca > 0:
            epi_result = cmd.super(aln_epi_sel, ref_epi_sel, cycles=0)
            epitope_rmsd = epi_result[0]
            epitope_atoms = epi_result[1]
        else:
            epitope_rmsd = -1.0
            epitope_atoms = 0
    except Exception as e:
        epitope_rmsd = -1.0
        epitope_atoms = 0

    print(f"{pdb_id:<8} {aln_epi_ca:<8} {global_align_rmsd:>12.3f}A {global_super_rmsd:>12.3f}A {super_atoms:>10} {epitope_rmsd:>10.3f}A {epitope_atoms:>8}")

    results[pdb_id] = {
        "global_align_rmsd": global_align_rmsd,
        "global_super_rmsd": global_super_rmsd,
        "super_atoms": super_atoms,
        "epitope_rmsd": epitope_rmsd,
        "epitope_atoms": epitope_atoms,
        "ref_epitope_ca": ref_epi_ca,
        "aligned_epitope_ca": aln_epi_ca,
    }

    # Keep structure for final visualization
    cmd.delete(obj_name)
    cmd.load(str(cif_path), obj_name)

print("="*90)

# Save results to JSON
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {output_file}")

# Analysis
print("\n" + "="*90)
print("ANALYSIS")
print("="*90)

print("""
Key observations:

1. These are DIFFERENT proteins grouped together based on ESM-2 embedding similarity.
   They share similar epitope "signatures" in embedding space, but NOT structural homology.

2. The "Epitope RMSD" calculated by cmd.super() uses structure-based superposition.
   - When epitopes are structurally similar (similar 3D shape), RMSD will be low.
   - When epitopes are NOT structurally similar, PyMOL matches few atoms → low RMSD but few atoms.

3. The small number of aligned epitope atoms (EpiAtoms) indicates that these proteins
   do NOT share similar epitope structures in 3D space.

4. For your goal of designing common antibodies:
   - ESM-2 embeddings capture SEQUENCE similarity/patterns, not 3D structure.
   - Low epitope RMSD with few atoms means the structures are NOT well aligned.
   - You need epitopes with both similar 3D structure AND orientation.

5. RECOMMENDATION:
   - Consider using structure-based alignment BEFORE grouping.
   - Or use RMSD threshold on epitope atoms as an additional grouping criterion.
   - Current grouping is based on embedding similarity alone.
""")

# Visualize
cmd.show("cartoon", "all")
cmd.color("gray80", "ref")

# Color aligned structures
colors = ["red", "orange", "yellow", "green", "cyan", "violet", "magenta", "salmon", "lime"]
for i, cif_path in enumerate(aligned_files):
    pdb_id = cif_path.stem.split('_')[0]
    obj_name = f"aln_{pdb_id}"
    cmd.color(colors[i % len(colors)], obj_name)

cmd.zoom("all")
cmd.bg_color("white")

# Save session
session_file = str(base_dir / "alignment_check_v2.pse")
cmd.save(session_file)
print(f"\nPyMOL session saved to {session_file}")
