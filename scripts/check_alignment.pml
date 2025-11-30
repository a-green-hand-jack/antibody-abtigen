# PyMOL script to check alignment quality
# Run with: mamba run -n pymol-env pymol -cq check_alignment.pml

import json
from pymol import cmd

# Paths (adjust as needed)
base_dir = "/ibex/user/wuj0c/Projects/Protein/GRADIENT/antibody-abtigen/data/epitope_pipeline/aligned_test/group_0000"
ref_cif = base_dir + "/reference/1A14_antigen.cif"
output_file = base_dir + "/alignment_rmsd_results.json"

# Aligned structures
aligned_files = [
    ("1A2Y", base_dir + "/aligned/1A2Y_antigen.cif"),
    ("1ADQ", base_dir + "/aligned/1ADQ_antigen.cif"),
    ("1AFV", base_dir + "/aligned/1AFV_antigen.cif"),
    ("1AHW", base_dir + "/aligned/1AHW_antigen.cif"),
    ("1AR1", base_dir + "/aligned/1AR1_antigen.cif"),
    ("1BGX", base_dir + "/aligned/1BGX_antigen.cif"),
    ("1BJ1", base_dir + "/aligned/1BJ1_antigen.cif"),
    ("1BQL", base_dir + "/aligned/1BQL_antigen.cif"),
    ("1BVK", base_dir + "/aligned/1BVK_antigen.cif"),
]

# Epitope residue selections (from epitope_residues.csv)
epitope_selections = {
    "1A14": "chain N and (resi 327 or resi 328 or resi 329 or resi 330 or resi 331 or resi 332 or resi 341 or resi 342 or resi 343 or resi 344 or resi 366 or resi 367 or resi 368 or resi 369 or resi 370 or resi 372 or resi 400 or resi 401 or resi 402 or resi 403 or resi 432)",
    "1A2Y": "chain C and (resi 18 or resi 19 or resi 22 or resi 23 or resi 24 or resi 25 or resi 26 or resi 27 or resi 102 or resi 103 or resi 116 or resi 117 or resi 118 or resi 119 or resi 120 or resi 121 or resi 122 or resi 124 or resi 125 or resi 128)",
    "1ADQ": "chain A and (resi 251 or resi 252 or resi 253 or resi 254 or resi 255 or resi 311 or resi 384 or resi 385 or resi 386 or resi 422 or resi 424 or resi 428 or resi 433 or resi 434 or resi 435 or resi 436 or resi 437 or resi 438 or resi 440)",
    "1AFV": "(chain A and (resi 71 or resi 72 or resi 74 or resi 75 or resi 76 or resi 77 or resi 78 or resi 79 or resi 81 or resi 82 or resi 83 or resi 84 or resi 85 or resi 100 or resi 101 or resi 102 or resi 103 or resi 136)) or (chain B and (resi 72 or resi 74 or resi 75 or resi 76 or resi 77 or resi 78 or resi 79 or resi 81 or resi 82 or resi 83 or resi 84 or resi 85 or resi 100 or resi 101 or resi 102 or resi 136))",
    "1AHW": "chain C and (resi 149 or resi 150 or resi 151 or resi 152 or resi 154)",
    "1AR1": "chain C and (resi 32 or resi 33 or resi 34 or resi 35 or resi 37 or resi 38 or resi 45 or resi 46 or resi 47 or resi 48 or resi 49 or resi 50 or resi 51 or resi 52 or resi 60 or resi 75 or resi 76 or resi 77 or resi 78 or resi 79 or resi 94 or resi 95)",
    "1BGX": "chain C and (resi 59 or resi 60 or resi 61 or resi 62 or resi 63 or resi 64 or resi 65 or resi 66 or resi 67 or resi 89 or resi 90 or resi 93 or resi 96 or resi 97 or resi 99 or resi 100 or resi 116 or resi 118 or resi 119 or resi 120 or resi 121 or resi 125)",
    "1BJ1": "chain P and (resi 113 or resi 114 or resi 116 or resi 117)",
    "1BQL": "chain Y and (resi 33 or resi 34 or resi 35 or resi 36 or resi 37 or resi 38 or resi 39 or resi 40 or resi 41 or resi 42)",
    "1BVK": "chain C and (resi 19 or resi 20 or resi 21 or resi 22 or resi 23 or resi 24 or resi 41 or resi 89 or resi 90 or resi 91 or resi 102 or resi 103 or resi 117 or resi 118 or resi 119 or resi 120 or resi 121)",
}

results = {}

# Load reference
cmd.load(ref_cif, "ref")
ref_epi_sel = epitope_selections["1A14"]
print(f"Reference: 1A14")
print(f"Reference epitope: {ref_epi_sel}")
print(f"Reference epitope CA atoms: {cmd.count_atoms('ref and (' + ref_epi_sel + ') and name CA')}")

print("\n" + "="*70)
print(f"{'PDB':<8} {'Global(Align)':<14} {'Global(Super)':<14} {'SuperAtoms':<12} {'Epitope':<12} {'EpiAtoms'}")
print("-"*70)

for pdb_id, cif_path in aligned_files:
    obj_name = f"aln_{pdb_id}"

    # Load aligned structure
    cmd.load(cif_path, obj_name)

    aln_epi_sel = epitope_selections.get(pdb_id, "none")

    # Calculate global RMSD using align (sequence-based, calculates RMSD without moving)
    try:
        # align returns (RMSD, aligned_atoms, cycles, RMSD_end, aligned_atoms_end, score, paired_residues)
        align_result = cmd.align(f"{obj_name} and name CA", "ref and name CA", cycles=0)
        global_align_rmsd = align_result[0]
    except Exception as e:
        print(f"Global align failed for {pdb_id}: {e}")
        global_align_rmsd = -1.0

    # Calculate global RMSD using super (structure-based, sequence-independent)
    # First reload to reset position
    cmd.delete(obj_name)
    cmd.load(cif_path, obj_name)

    try:
        super_result = cmd.super(f"{obj_name} and name CA", "ref and name CA", cycles=0)
        global_super_rmsd = super_result[0]
        super_atoms = super_result[1]
    except Exception as e:
        print(f"Super failed for {pdb_id}: {e}")
        global_super_rmsd = -1.0
        super_atoms = 0

    # Calculate epitope-only RMSD
    # Reload again to reset position
    cmd.delete(obj_name)
    cmd.load(cif_path, obj_name)

    try:
        epi_sel_aln = f"{obj_name} and ({aln_epi_sel}) and name CA"
        epi_sel_ref = f"ref and ({ref_epi_sel}) and name CA"

        aln_epi_atoms = cmd.count_atoms(epi_sel_aln)
        ref_epi_atoms = cmd.count_atoms(epi_sel_ref)

        if aln_epi_atoms > 0 and ref_epi_atoms > 0:
            epi_result = cmd.super(epi_sel_aln, epi_sel_ref, cycles=0)
            epitope_rmsd = epi_result[0]
            epitope_atoms = epi_result[1]
        else:
            epitope_rmsd = -1.0
            epitope_atoms = 0
    except Exception as e:
        print(f"Epitope super failed for {pdb_id}: {e}")
        epitope_rmsd = -1.0
        epitope_atoms = 0

    print(f"{pdb_id:<8} {global_align_rmsd:>12.3f}A {global_super_rmsd:>12.3f}A {super_atoms:>10} {epitope_rmsd:>10.3f}A {epitope_atoms:>8}")

    results[pdb_id] = {
        "global_align_rmsd": global_align_rmsd,
        "global_super_rmsd": global_super_rmsd,
        "super_atoms": super_atoms,
        "epitope_rmsd": epitope_rmsd,
        "epitope_atoms": epitope_atoms,
    }

    # Keep the structure for visualization (reload once more)
    cmd.delete(obj_name)
    cmd.load(cif_path, obj_name)

print("="*70)

# Save results to JSON
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {output_file}")

# Visualize
cmd.show("cartoon", "all")
cmd.color("gray80", "ref")

# Highlight epitope on reference
cmd.color("marine", f"ref and ({ref_epi_sel})")
cmd.show("sticks", f"ref and ({ref_epi_sel})")

# Color aligned structures
colors = ["red", "orange", "yellow", "green", "cyan", "violet", "magenta", "salmon", "lime"]
for i, (pdb_id, _) in enumerate(aligned_files):
    obj_name = f"aln_{pdb_id}"
    cmd.color(colors[i % len(colors)], obj_name)
    # Highlight epitope
    epi_sel = epitope_selections.get(pdb_id, "none")
    if epi_sel != "none":
        cmd.show("sticks", f"{obj_name} and ({epi_sel})")

cmd.zoom("all")
cmd.bg_color("white")

# Save session
session_file = base_dir + "/alignment_check.pse"
cmd.save(session_file)
print(f"\nPyMOL session saved to {session_file}")
