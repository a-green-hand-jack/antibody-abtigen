import argparse
import json
import logging
import shutil
from pathlib import Path
from typing import List, Set

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Copy test YAML files to separate test directories.")
    parser.add_argument("--test_entries", type=str, default="data/split/test_entry.json", help="Path to test entry JSON file.")
    parser.add_argument("--single_antigen_dir", type=str, default="data/yamls/boltzgen_like/single_antigen", help="Source directory for single-antigen YAMLs.")
    parser.add_argument("--multi_antigen_dir", type=str, default="data/yamls/boltzgen_like/multi_antigen", help="Source directory for multi-antigen YAMLs.")
    parser.add_argument("--single_antigen_test_dir", type=str, default="data/yamls/boltzgen_like/single_antigen_tests", help="Output directory for single-antigen test YAMLs.")
    parser.add_argument("--multi_antigen_test_dir", type=str, default="data/yamls/boltzgen_like/multi_antigen_tests", help="Output directory for multi-antigen test YAMLs.")
    return parser.parse_args()


def setup_directories(dirs: List[Path]):
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def load_test_entries(test_entries_path: Path) -> Set[str]:
    """Load test entry identifiers from JSON file."""
    try:
        with open(test_entries_path, 'r') as f:
            entries = json.load(f)
        return set(entries)
    except Exception as e:
        logger.error(f"Failed to load test entries from {test_entries_path}: {e}")
        return set()


def extract_pdb_chain_from_yaml_name(yaml_filename: str) -> str:
    """
    Extract PDB_chain from YAML filename.
    For single_antigen: 8a1e_A.yaml -> 8a1e_A
    For multi_antigen: cluster_1.yaml -> need to match by content
    """
    # Remove .yaml extension
    return yaml_filename.replace('.yaml', '')


def find_matching_pdbs_in_multi_antigen(yaml_content: str, test_entries: Set[str]) -> bool:
    """
    Check if any PDB entry from test set appears in multi-antigen YAML content.
    Multi-antigen YAMLs contain file paths like: ../../../antibody_antigen_single/8a1e_A_antibody.cif
    """
    for entry in test_entries:
        pdb_chain = entry.split('_')[0:3]
        pdb_code = pdb_chain[0]  # e.g., "8a1e" from "8a1e_H_L_A"

        if f"{pdb_code}_" in yaml_content:
            return True

    return False


def main():
    args = parse_arguments()

    single_antigen_dir = Path(args.single_antigen_dir)
    multi_antigen_dir = Path(args.multi_antigen_dir)
    single_antigen_test_dir = Path(args.single_antigen_test_dir)
    multi_antigen_test_dir = Path(args.multi_antigen_test_dir)

    setup_directories([single_antigen_test_dir, multi_antigen_test_dir])

    # Load test entries
    logger.info(f"Loading test entries from {args.test_entries}")
    test_entries = load_test_entries(Path(args.test_entries))
    logger.info(f"Found {len(test_entries)} test entries")

    # Extract PDB codes from test entries (e.g., "8a1e" from "8a1e_H_L_A")
    test_pdb_codes = set()
    for entry in test_entries:
        pdb_code = entry.split('_')[0]
        test_pdb_codes.add(pdb_code)

    logger.info(f"Test PDB codes: {len(test_pdb_codes)}")

    # Process single-antigen YAMLs
    logger.info(f"Processing single-antigen YAMLs from {single_antigen_dir}")
    single_antigen_count = 0
    if single_antigen_dir.exists():
        for yaml_file in sorted(single_antigen_dir.glob('*.yaml')):
            # Extract PDB_chain from filename (e.g., "8a1e_A" from "8a1e_A.yaml")
            pdb_chain_str = extract_pdb_chain_from_yaml_name(yaml_file.name)
            pdb_code = pdb_chain_str.split('_')[0]

            # Check if this PDB code is in test set
            if pdb_code in test_pdb_codes:
                # Copy YAML file
                dst_yaml = single_antigen_test_dir / yaml_file.name
                shutil.copy2(yaml_file, dst_yaml)
                logger.info(f"Copied {yaml_file.name} to single-antigen tests")
                single_antigen_count += 1

    # Process multi-antigen YAMLs
    logger.info(f"Processing multi-antigen YAMLs from {multi_antigen_dir}")
    multi_antigen_count = 0
    if multi_antigen_dir.exists():
        for yaml_file in sorted(multi_antigen_dir.glob('*.yaml')):
            # Read the YAML file to check if it contains test PDbs
            try:
                with open(yaml_file, 'r') as f:
                    yaml_content = f.read()

                # Check if any test PDB code appears in the YAML content
                is_test_yaml = False
                for pdb_code in test_pdb_codes:
                    if f"{pdb_code}_" in yaml_content:
                        is_test_yaml = True
                        break

                if is_test_yaml:
                    # Copy YAML file
                    dst_yaml = multi_antigen_test_dir / yaml_file.name
                    shutil.copy2(yaml_file, dst_yaml)
                    logger.info(f"Copied {yaml_file.name} to multi-antigen tests")
                    multi_antigen_count += 1

            except Exception as e:
                logger.warning(f"Failed to process {yaml_file.name}: {e}")

    logger.info(f"\nSummary:")
    logger.info(f"  Single-antigen YAMLs copied: {single_antigen_count} -> {single_antigen_test_dir}")
    logger.info(f"  Multi-antigen YAMLs copied: {multi_antigen_count} -> {multi_antigen_test_dir}")


if __name__ == "__main__":
    main()
