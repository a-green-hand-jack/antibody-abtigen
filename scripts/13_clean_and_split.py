import argparse
import logging
from pathlib import Path
from typing import Set

import gemmi
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Clean and split aligned structures into antigen and antibody files."
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        required=True,
        help="Path to the epitope cluster summary CSV file.",
    )
    parser.add_argument(
        "--input_aligned_dir",
        type=str,
        required=True,
        help="Directory containing the aligned CIF files (e.g., data/aligned).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Base directory to save cleaned and split files (e.g., data/cleaned_split).",
    )
    parser.add_argument(
        "--workers", type=int, default=8, help="Number of parallel workers."
    )
    return parser.parse_args()


def parse_antigen_chains(ag_chain_str: str) -> Set[str]:
    """
    Parses the antigen chain string (e.g., "A|B", "['A', 'B']", "A") into a set of chain IDs.
    """
    if pd.isna(ag_chain_str) or not ag_chain_str or ag_chain_str == "nan":
        return set()

    clean_str = (
        str(ag_chain_str)
        .replace("|", ",")
        .replace("[", "")
        .replace("]", "")
        .replace("'", "")
        .replace('"', "")
    )
    parts = [p.strip() for p in clean_str.split(",") if p.strip()]
    return set(parts)


def process_entry(row: pd.Series, input_aligned_dir: Path, output_base_dir: Path):
    """
    Cleans and splits a single aligned structure entry.
    """
    cluster_id = row["cluster_id"]
    epitope_id = row["epitope_id"]

    # Skip unclustered if necessary, but usually we want to process everything in the summary
    if cluster_id == -1:
        return

    # Path to the ALIGNED input file
    input_path = input_aligned_dir / str(cluster_id) / f"{epitope_id}.cif"

    if not input_path.exists():
        # This might happen if the cluster size was < min_size in the previous step
        # logger.debug(f"Input file not found (likely skipped in alignment): {input_path}")
        return

    try:
        # 1. Read Structure
        st = gemmi.read_structure(str(input_path))

        # 2. Clean: Remove Waters
        st.remove_waters()
        # st.remove_ligands_and_waters() # Maybe too aggressive if ligand is part of antigen?
        # Let's stick to removing waters and then filtering by chain ID.

        # 3. Identify Target Chains
        h_chain = str(row["H_chain"]) if pd.notna(row["H_chain"]) else ""
        l_chain = str(row["L_chain"]) if pd.notna(row["L_chain"]) else ""
        ag_chains = parse_antigen_chains(row["antigen_chains"])

        ab_chains = {h_chain, l_chain} - {""}

        # 4. Create New Structures
        st_ab = gemmi.Structure()
        st_ab.name = f"{epitope_id}_antibody"
        model_ab = gemmi.Model("1")
        st_ab.add_model(model_ab)

        st_ag = gemmi.Structure()
        st_ag.name = f"{epitope_id}_antigen"
        model_ag = gemmi.Model("1")
        st_ag.add_model(model_ag)

        # 5. Distribute Chains
        # We iterate over the first model (assuming standard PDB structure)
        source_model = st[0]

        found_ab = False
        found_ag = False

        for chain in source_model:
            # Check against chain.name (which is label_asym_id in mmCIF usually, or auth_asym_id depending on how Gemmi read it)
            # Gemmi read_structure usually maps label_asym_id to chain.name by default for CIF.
            # Our metadata usually comes from SAbDab summary which uses auth_asym_id?
            # This is a potential risk point.
            # However, 09_extract_epitope used `model[id]` which looks up by name.
            # Let's try to match both if possible, or trust `chain.name`.

            cid = chain.name

            if cid in ab_chains:
                model_ab.add_chain(chain.clone())
                found_ab = True
            elif cid in ag_chains:
                model_ag.add_chain(chain.clone())
                found_ag = True
            # Else: Ignore (non-interacting chain or junk)

        # 6. Save
        if found_ab and found_ag:
            output_dir = output_base_dir / str(cluster_id) / epitope_id
            output_dir.mkdir(parents=True, exist_ok=True)

            st_ab.make_mmcif_document().write_file(str(output_dir / "antibody.cif"))
            st_ag.make_mmcif_document().write_file(str(output_dir / "antigen.cif"))
        else:
            logger.warning(
                f"Skipping {epitope_id}: Missing chains (Found Ab: {found_ab}, Found Ag: {found_ag})"
            )

    except Exception as e:
        logger.error(f"Error processing {epitope_id}: {e}")


def main():
    args = parse_arguments()

    input_csv_path = Path(args.input_csv)
    input_aligned_dir = Path(args.input_aligned_dir)
    output_dir = Path(args.output_dir)

    if not input_csv_path.exists():
        logger.error(f"Input CSV not found: {input_csv_path}")
        return

    df = pd.read_csv(input_csv_path)
    logger.info(f"Loaded {len(df)} records.")

    # Parallel Processing
    Parallel(n_jobs=args.workers)(
        delayed(process_entry)(row, input_aligned_dir, output_dir)
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Cleaning and Splitting")
    )

    logger.info("Cleaning and splitting complete.")


if __name__ == "__main__":
    main()
