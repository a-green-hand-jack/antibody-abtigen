import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import gemmi
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Water residue names to exclude
WATER_RESIDUES = {"HOH", "WAT", "DOD", "H2O"}

# Categories that need filtering by chain (auth_asym_id)
CHAIN_FILTER_CATEGORIES = {
    "_atom_site": "auth_asym_id",
    "_pdbx_poly_seq_scheme": "pdb_strand_id",
    "_pdbx_nonpoly_scheme": "pdb_strand_id",
    "_struct_site_gen": "auth_asym_id",
}

# Categories that need filtering by entity_id
ENTITY_FILTER_CATEGORIES = {
    "_entity": "id",
    "_entity_poly": "entity_id",
    "_entity_poly_seq": "entity_id",
    "_entity_src_gen": "entity_id",
    "_entity_src_nat": "entity_id",
    "_pdbx_entity_nonpoly": "entity_id",
    "_pdbx_entity_branch": "entity_id",
    "_pdbx_entity_branch_list": "entity_id",
    "_pdbx_entity_branch_link": "entity_id",
}

# Categories that need special handling (multiple chain columns)
MULTI_CHAIN_CATEGORIES = {
    "_struct_conn": ["ptnr1_auth_asym_id", "ptnr2_auth_asym_id"],
}

# Categories that need filtering by label_asym_id (internal chain ID)
LABEL_ASYM_FILTER_CATEGORIES = {
    "_struct_asym": "id",
    "_pdbx_branch_scheme": "asym_id",
    "_pdbx_nonpoly_scheme": "asym_id",
}


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


def get_chain_entity_mapping(block: gemmi.cif.Block) -> Dict[str, str]:
    """
    从 _atom_site 获取 auth_asym_id -> label_entity_id 的映射。
    Returns: {"A": "1", "B": "2", "C": "5", ...}
    """
    mapping = {}
    try:
        atom_site = block.find_mmcif_category("_atom_site.")
        if not atom_site:
            return mapping

        chain_col = atom_site.find_column("auth_asym_id")
        entity_col = atom_site.find_column("label_entity_id")

        for chain_id, entity_id in zip(chain_col, entity_col):
            chain_id = chain_id.strip()
            entity_id = entity_id.strip()
            if chain_id and chain_id not in mapping:
                mapping[chain_id] = entity_id
    except (ValueError, RuntimeError):
        pass

    return mapping


def get_label_asym_id_mapping(block: gemmi.cif.Block) -> Dict[str, Tuple[str, str]]:
    """
    从 _atom_site 获取 label_asym_id -> (auth_asym_id, label_entity_id) 的映射。
    Returns: {"A": ("A", "1"), "EA": ("A", "11"), ...}

    label_asym_id 是内部链 ID，每个结构单元唯一。
    auth_asym_id 是作者链 ID，蛋白和其附属配体共用。
    """
    mapping = {}
    try:
        atom_site = block.find_mmcif_category("_atom_site.")
        if not atom_site:
            return mapping

        label_asym_col = atom_site.find_column("label_asym_id")
        auth_asym_col = atom_site.find_column("auth_asym_id")
        entity_col = atom_site.find_column("label_entity_id")

        for label_asym, auth_asym, entity_id in zip(label_asym_col, auth_asym_col, entity_col):
            label_asym = label_asym.strip()
            auth_asym = auth_asym.strip()
            entity_id = entity_id.strip()
            if label_asym and label_asym not in mapping:
                mapping[label_asym] = (auth_asym, entity_id)
    except (ValueError, RuntimeError):
        pass

    return mapping


def get_connected_ligand_info(
    block: gemmi.cif.Block, target_auth_chain_ids: Set[str]
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    从 _struct_conn 找出与目标链共价连接的配体信息。

    返回: (connected_auth_asym_ids, connected_label_asym_ids, connected_entity_ids)

    在 mmCIF 中，蛋白链 A 可能有附属的糖链 (如 NAG)，它们的:
    - auth_asym_id 可能相同 (都是 "A") 或不同 (如 "K")
    - label_asym_id 不同 (蛋白是 "A"，糖是 "EA" 或 "S")
    - entity_id 不同 (蛋白是 "1"，糖是 "11" 或 "7")

    我们需要找出这些共价连接的配体，保留它们的原子和 entity 定义。
    """
    connected_auth_asym_ids: Set[str] = set()
    connected_label_asym_ids: Set[str] = set()
    connected_entity_ids: Set[str] = set()

    try:
        # 获取 label_asym_id 映射
        label_asym_mapping = get_label_asym_id_mapping(block)

        # 找出目标链对应的 label_asym_ids
        target_label_asym_ids = set()
        for label_asym, (auth_asym, entity_id) in label_asym_mapping.items():
            if auth_asym in target_auth_chain_ids:
                target_label_asym_ids.add(label_asym)

        # 查找 _struct_conn 中的共价连接
        struct_conn = block.find_mmcif_category("_struct_conn.")
        if not struct_conn:
            return connected_auth_asym_ids, connected_label_asym_ids, connected_entity_ids

        # 获取列索引
        try:
            conn_type_col = struct_conn.find_column("conn_type_id")
            ptnr1_label_asym_col = struct_conn.find_column("ptnr1_label_asym_id")
            ptnr2_label_asym_col = struct_conn.find_column("ptnr2_label_asym_id")
        except (ValueError, RuntimeError):
            return connected_auth_asym_ids, connected_label_asym_ids, connected_entity_ids

        # 遍历所有连接
        for conn_type, ptnr1_asym, ptnr2_asym in zip(
            conn_type_col, ptnr1_label_asym_col, ptnr2_label_asym_col
        ):
            conn_type = conn_type.strip()
            ptnr1_asym = ptnr1_asym.strip()
            ptnr2_asym = ptnr2_asym.strip()

            # 只处理共价连接 (covale, covale_base, covale_phosph, covale_sugar)
            if not conn_type.lower().startswith("covale"):
                continue

            # 如果一端是目标链，添加另一端
            if ptnr1_asym in target_label_asym_ids and ptnr2_asym not in target_label_asym_ids:
                connected_label_asym_ids.add(ptnr2_asym)
                if ptnr2_asym in label_asym_mapping:
                    auth_asym, entity_id = label_asym_mapping[ptnr2_asym]
                    connected_auth_asym_ids.add(auth_asym)
                    connected_entity_ids.add(entity_id)

            if ptnr2_asym in target_label_asym_ids and ptnr1_asym not in target_label_asym_ids:
                connected_label_asym_ids.add(ptnr1_asym)
                if ptnr1_asym in label_asym_mapping:
                    auth_asym, entity_id = label_asym_mapping[ptnr1_asym]
                    connected_auth_asym_ids.add(auth_asym)
                    connected_entity_ids.add(entity_id)

    except (ValueError, RuntimeError):
        pass

    return connected_auth_asym_ids, connected_label_asym_ids, connected_entity_ids


def get_entity_ids_for_chains(
    block: gemmi.cif.Block, chain_ids: Set[str]
) -> Set[str]:
    """
    根据 auth_asym_id 找到对应的所有 entity_id 集合。

    注意: 一个 auth_asym_id 可能对应多个 entity_id，
    例如蛋白链 A (entity=1) 和其附属的 NAG 配体 (entity=11) 都有 auth_asym_id=A。
    """
    entity_ids = set()
    try:
        atom_site = block.find_mmcif_category("_atom_site.")
        if not atom_site:
            return entity_ids

        chain_col = atom_site.find_column("auth_asym_id")
        entity_col = atom_site.find_column("label_entity_id")

        for chain_id, entity_id in zip(chain_col, entity_col):
            chain_id = chain_id.strip()
            entity_id = entity_id.strip()
            if chain_id in chain_ids and entity_id:
                entity_ids.add(entity_id)
    except (ValueError, RuntimeError):
        pass

    return entity_ids


def get_category_from_tag(tag: str) -> str:
    """从 tag 中提取 category 名称。例如 '_atom_site.id' -> '_atom_site'"""
    if "." in tag:
        return tag.rsplit(".", 1)[0]
    return tag


def get_loop_rows(loop: gemmi.cif.Loop) -> List[List[str]]:
    """从 gemmi.cif.Loop 对象中提取所有行。"""
    rows = []
    n_cols = loop.width()
    n_rows = loop.length()
    for i in range(n_rows):
        row = [loop[i, j] for j in range(n_cols)]
        rows.append(row)
    return rows


def filter_loop_rows(
    loop: gemmi.cif.Loop,
    filter_col_name: str,
    allowed_values: Set[str],
    exclude_water: bool = False,
    comp_id_col_name: Optional[str] = None,
) -> Tuple[List[str], List[List[str]]]:
    """
    过滤 loop 的行，返回 (tags, filtered_rows)。

    Args:
        loop: gemmi.cif.Loop 对象
        filter_col_name: 用于过滤的列名（不含 category 前缀）
        allowed_values: 允许的值集合
        exclude_water: 是否排除水分子
        comp_id_col_name: 残基类型列名（用于排除水）
    """
    tags = list(loop.tags)

    # 找到过滤列的索引
    filter_idx = None
    comp_id_idx = None

    for i, tag in enumerate(tags):
        col_name = tag.split(".")[-1] if "." in tag else tag
        if col_name == filter_col_name:
            filter_idx = i
        if comp_id_col_name and col_name == comp_id_col_name:
            comp_id_idx = i

    # 获取所有行
    all_rows = get_loop_rows(loop)

    if filter_idx is None:
        # 找不到过滤列，返回所有行
        return tags, all_rows

    filtered_rows = []
    for row_list in all_rows:
        val = row_list[filter_idx].strip()

        # 检查是否在允许的值中
        if val not in allowed_values:
            continue

        # 检查是否需要排除水
        if exclude_water and comp_id_idx is not None:
            comp_id = row_list[comp_id_idx].strip()
            if comp_id in WATER_RESIDUES:
                continue

        filtered_rows.append(row_list)

    return tags, filtered_rows


def filter_loop_multi_chain(
    loop: gemmi.cif.Loop,
    filter_col_names: List[str],
    allowed_values: Set[str],
) -> Tuple[List[str], List[List[str]]]:
    """
    过滤 loop 的行，要求所有指定列的值都在 allowed_values 中。
    用于 _struct_conn 等有多个链列的 category。
    """
    tags = list(loop.tags)

    # 找到所有过滤列的索引
    filter_indices = []
    for i, tag in enumerate(tags):
        col_name = tag.split(".")[-1] if "." in tag else tag
        if col_name in filter_col_names:
            filter_indices.append(i)

    # 获取所有行
    all_rows = get_loop_rows(loop)

    if not filter_indices:
        return tags, all_rows

    filtered_rows = []
    for row_list in all_rows:
        # 检查所有过滤列的值是否都在允许集合中
        all_match = True
        for idx in filter_indices:
            val = row_list[idx].strip()
            if val and val not in (".", "?") and val not in allowed_values:
                all_match = False
                break
        if all_match:
            filtered_rows.append(row_list)

    return tags, filtered_rows


def create_filtered_block(
    source_block: gemmi.cif.Block,
    target_chain_ids: Set[str],
    target_entity_ids: Set[str],
    target_label_asym_ids: Set[str],
    new_entry_name: str,
) -> gemmi.cif.Block:
    """
    复制 source_block，过滤只保留 target_chain_ids 相关的数据。

    Args:
        source_block: 源 CIF block
        target_chain_ids: 要保留的链 ID 集合 (auth_asym_id)
        target_entity_ids: 要保留的 entity ID 集合
        target_label_asym_ids: 要保留的 label_asym_id 集合 (内部链 ID)
        new_entry_name: 新 block 的名称
    """
    # 创建新的 Document 和 Block
    new_block = gemmi.cif.Block(new_entry_name)

    # 遍历源 block 的所有 items
    for item in source_block:
        if item.loop is not None:
            # 这是一个 loop
            loop = item.loop
            if not loop.tags:
                continue

            category = get_category_from_tag(loop.tags[0])

            # 决定如何处理这个 loop
            if category in CHAIN_FILTER_CATEGORIES:
                filter_col = CHAIN_FILTER_CATEGORIES[category]
                exclude_water = category == "_atom_site"
                comp_id_col = "label_comp_id" if exclude_water else None

                tags, rows = filter_loop_rows(
                    loop,
                    filter_col,
                    target_chain_ids,
                    exclude_water=exclude_water,
                    comp_id_col_name=comp_id_col,
                )

                if rows:
                    # 创建新 loop
                    new_loop = new_block.init_loop("", tags)
                    for row in rows:
                        new_loop.add_row(row)

            elif category in ENTITY_FILTER_CATEGORIES:
                filter_col = ENTITY_FILTER_CATEGORIES[category]
                tags, rows = filter_loop_rows(
                    loop, filter_col, target_entity_ids
                )

                if rows:
                    new_loop = new_block.init_loop("", tags)
                    for row in rows:
                        new_loop.add_row(row)

            elif category in LABEL_ASYM_FILTER_CATEGORIES:
                # 使用 label_asym_id 过滤 (如 _struct_asym)
                filter_col = LABEL_ASYM_FILTER_CATEGORIES[category]
                tags, rows = filter_loop_rows(
                    loop, filter_col, target_label_asym_ids
                )

                if rows:
                    new_loop = new_block.init_loop("", tags)
                    for row in rows:
                        new_loop.add_row(row)

            elif category in MULTI_CHAIN_CATEGORIES:
                filter_cols = MULTI_CHAIN_CATEGORIES[category]
                tags, rows = filter_loop_multi_chain(
                    loop, filter_cols, target_chain_ids
                )

                if rows:
                    new_loop = new_block.init_loop("", tags)
                    for row in rows:
                        new_loop.add_row(row)

            else:
                # 直接复制整个 loop
                all_rows = get_loop_rows(loop)
                if all_rows:
                    new_loop = new_block.init_loop("", list(loop.tags))
                    for row in all_rows:
                        new_loop.add_row(row)

        elif item.pair is not None:
            # 这是一个单值 item (key-value pair)
            tag, value = item.pair
            # 更新 entry.id
            if tag == "_entry.id":
                new_block.set_pair(tag, gemmi.cif.quote(new_entry_name))
            else:
                new_block.set_pair(tag, value)

    return new_block


def process_entry(row: pd.Series, input_aligned_dir: Path, output_base_dir: Path):
    """
    Cleans and splits a single aligned structure entry.
    Preserves all mmCIF metadata from the original file.
    """
    cluster_id = row["cluster_id"]
    epitope_id = row["epitope_id"]

    # Skip unclustered if necessary
    if cluster_id == -1:
        return

    # Path to the ALIGNED input file
    input_path = input_aligned_dir / str(cluster_id) / f"{epitope_id}.cif"

    if not input_path.exists():
        return

    try:
        # 1. Read CIF as Document to preserve all metadata
        doc = gemmi.cif.read_file(str(input_path))
        block = doc.sole_block()

        # 2. Identify Target Chains
        h_chain = str(row["H_chain"]) if pd.notna(row["H_chain"]) else ""
        l_chain = str(row["L_chain"]) if pd.notna(row["L_chain"]) else ""
        ag_chains = parse_antigen_chains(row["antigen_chains"])

        ab_chains = {h_chain, l_chain} - {""}

        # 3. Validate chains exist in the structure
        chain_entity_map = get_chain_entity_mapping(block)
        existing_chains = set(chain_entity_map.keys())

        ab_chains = ab_chains & existing_chains
        ag_chains = ag_chains & existing_chains

        if not ab_chains or not ag_chains:
            return

        # 4. Get entity IDs for each chain set
        ab_entity_ids = get_entity_ids_for_chains(block, ab_chains)
        ag_entity_ids = get_entity_ids_for_chains(block, ag_chains)

        # 5. Get label_asym_ids and find connected ligands
        label_asym_mapping = get_label_asym_id_mapping(block)

        # 5a. 找出 ab_chains 对应的 label_asym_ids
        ab_label_asym_ids = set()
        for label_asym, (auth_asym, entity_id) in label_asym_mapping.items():
            if auth_asym in ab_chains:
                ab_label_asym_ids.add(label_asym)

        # 5b. 找出 ag_chains 对应的 label_asym_ids
        ag_label_asym_ids = set()
        for label_asym, (auth_asym, entity_id) in label_asym_mapping.items():
            if auth_asym in ag_chains:
                ag_label_asym_ids.add(label_asym)

        # 5c. 找出共价连接的配体，扩展 chain_ids, entity_ids 和 label_asym_ids
        ab_conn_auth, ab_conn_label, ab_conn_entity = get_connected_ligand_info(
            block, ab_chains
        )
        ag_conn_auth, ag_conn_label, ag_conn_entity = get_connected_ligand_info(
            block, ag_chains
        )

        # 扩展 auth_asym_ids (用于 _atom_site 等按 auth_asym_id 过滤的 category)
        ab_chains = ab_chains | ab_conn_auth
        ag_chains = ag_chains | ag_conn_auth

        ab_entity_ids = ab_entity_ids | ab_conn_entity
        ag_entity_ids = ag_entity_ids | ag_conn_entity
        ab_label_asym_ids = ab_label_asym_ids | ab_conn_label
        ag_label_asym_ids = ag_label_asym_ids | ag_conn_label

        # 6. Create filtered blocks
        ab_block = create_filtered_block(
            block, ab_chains, ab_entity_ids, ab_label_asym_ids, f"{epitope_id}_antibody"
        )
        ag_block = create_filtered_block(
            block, ag_chains, ag_entity_ids, ag_label_asym_ids, f"{epitope_id}_antigen"
        )

        # 7. Save
        output_dir = output_base_dir / str(cluster_id) / epitope_id
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write antibody CIF
        ab_doc = gemmi.cif.Document()
        ab_doc.add_copied_block(ab_block)
        ab_doc.write_file(str(output_dir / "antibody.cif"))

        # Write antigen CIF
        ag_doc = gemmi.cif.Document()
        ag_doc.add_copied_block(ag_block)
        ag_doc.write_file(str(output_dir / "antigen.cif"))

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
