"""三元组对齐输出模块 - 使用 PyMOL 进行对齐并保存完整 CIF"""
import os
import subprocess
import json
import pandas as pd
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def _get_pymol_python() -> Optional[str]:
    """获取 PyMOL Python 路径"""
    # macOS PyMOL.app
    pymol_app = "/Applications/PyMOL.app/Contents/bin/python"
    if os.path.exists(pymol_app):
        return pymol_app
    return None


def _align_and_save_triplet_pymol(
    human_ag_file: str,
    mouse_ag_file: str,
    antibody_file: str,
    output_dir: str,
    human_chain: str,
    mouse_chain: str,
) -> Tuple[float, bool]:
    """
    使用 PyMOL 对齐三元组并保存为完整 CIF 文件

    Args:
        human_ag_file: 人源抗原 CIF 文件
        mouse_ag_file: 鼠源抗原 CIF 文件（未对齐的原始文件）
        antibody_file: 抗体 CIF 文件
        output_dir: 输出目录
        human_chain: 人源抗原链 ID（用于对齐）
        mouse_chain: 鼠源抗原链 ID（用于对齐）

    Returns:
        (RMSD, success)
    """
    python_path = _get_pymol_python()
    if not python_path:
        raise RuntimeError("PyMOL not found. Please install PyMOL.app")

    os.makedirs(output_dir, exist_ok=True)

    human_out = os.path.join(output_dir, "human_antigen.cif")
    mouse_out = os.path.join(output_dir, "mouse_antigen.cif")
    antibody_out = os.path.join(output_dir, "antibody.cif")

    # PyMOL 脚本：加载、对齐、保存
    script = f'''
import sys
import json
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-qc'])

# 加载三个结构
cmd.load("{human_ag_file}", "human_ag")
cmd.load("{mouse_ag_file}", "mouse_ag")
cmd.load("{antibody_file}", "antibody")

# 对齐鼠源抗原到人源抗原
mouse_sel = "mouse_ag and chain {mouse_chain}"
human_sel = "human_ag and chain {human_chain}"

result = cmd.align(mouse_sel, human_sel)
rmsd = result[0] if isinstance(result, tuple) else result

# 保存为完整的 CIF 格式
cmd.save("{human_out}", "human_ag", format="cif")
cmd.save("{mouse_out}", "mouse_ag", format="cif")  # 对齐后的坐标
cmd.save("{antibody_out}", "antibody", format="cif")

# 输出结果
sys.stderr.write(json.dumps({{"rmsd": rmsd, "success": True}}) + "\\n")
sys.stderr.flush()
'''

    result = subprocess.run(
        [python_path, '-c', script],
        capture_output=True,
        text=True,
        timeout=120
    )

    if result.returncode != 0:
        logger.error(f"PyMOL failed: {result.stderr}")
        return 0.0, False

    # 解析结果
    for line in result.stderr.strip().split('\n'):
        line = line.strip()
        if line.startswith('{') and 'rmsd' in line:
            try:
                data = json.loads(line)
                return data['rmsd'], data['success']
            except json.JSONDecodeError:
                continue

    # 检查输出文件是否存在作为成功的备用判断
    if all(os.path.exists(f) for f in [human_out, mouse_out, antibody_out]):
        logger.warning("Could not parse RMSD but output files exist")
        return 0.0, True

    return 0.0, False


def generate_triplets(
    data_dir: str,
    output_dir: Optional[str] = None,
    dry_run: bool = False
) -> int:
    """
    从现有数据生成三元组对齐输出

    读取 human_mouse_pairs.csv，对每个配对：
    1. 加载人源抗原、鼠源抗原、抗体
    2. 使用 PyMOL 将鼠源抗原对齐到人源抗原
    3. 保存三个完整的 CIF 文件

    Args:
        data_dir: 数据目录 (包含 SAbDab/, MouseAntigen/, human_mouse_pairs.csv)
        output_dir: 输出目录，默认为 data_dir/HumanMouseAntigenAntibody
        dry_run: 只打印将要执行的操作，不实际执行

    Returns:
        成功生成的三元组数量
    """
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    if output_dir is None:
        output_dir = os.path.join(data_dir, "HumanMouseAntigenAntibody")

    # 读取映射表
    pairs_csv = os.path.join(data_dir, "human_mouse_pairs.csv")
    if not os.path.exists(pairs_csv):
        raise FileNotFoundError(f"Pairs CSV not found: {pairs_csv}")

    df = pd.read_csv(pairs_csv)
    success_count = 0
    total = len(df)

    logger.info(f"Processing {total} pairs from {pairs_csv}")
    logger.info(f"Output directory: {output_dir}")

    for idx, row in df.iterrows():
        triplet_name = f"{row['human_antigen_gene']}_{row['mouse_antigen_gene']}_{row['human_antigen_pdb']}"
        triplet_dir = os.path.join(output_dir, triplet_name)

        # 源文件路径
        human_ag_src = os.path.join(data_dir, row['human_antigen_file'])
        antibody_src = os.path.join(data_dir, row['antibody_file'])
        mouse_ag_src = os.path.join(data_dir, row['mouse_antigen_file'])  # 使用原始未对齐的鼠源抗原

        # 检查源文件是否存在
        missing = []
        for f, name in [(human_ag_src, 'human_ag'), (antibody_src, 'antibody'), (mouse_ag_src, 'mouse_ag')]:
            if not os.path.exists(f):
                missing.append(f"{name}: {f}")

        if missing:
            logger.warning(f"[{idx + 1}/{total}] Skipping {triplet_name}: missing files:")
            for m in missing:
                logger.warning(f"  - {m}")
            continue

        if dry_run:
            logger.info(f"[{idx + 1}/{total}] [DRY-RUN] Would create: {triplet_dir}")
            success_count += 1
            continue

        # 获取链 ID
        human_chains = row.get('human_antigen_chains', 'A')
        if isinstance(human_chains, str) and ',' in human_chains:
            human_chain = human_chains.split(',')[0].strip()
        else:
            human_chain = str(human_chains).strip() if pd.notna(human_chains) else 'A'

        mouse_chain = str(row.get('mouse_antigen_chain', 'A')).strip()

        # 使用 PyMOL 对齐并保存
        try:
            rmsd, ok = _align_and_save_triplet_pymol(
                human_ag_file=human_ag_src,
                mouse_ag_file=mouse_ag_src,
                antibody_file=antibody_src,
                output_dir=triplet_dir,
                human_chain=human_chain,
                mouse_chain=mouse_chain,
            )

            if ok:
                logger.info(f"[{idx + 1}/{total}] Created triplet: {triplet_name} (RMSD={rmsd:.2f}Å)")
                success_count += 1
            else:
                logger.error(f"[{idx + 1}/{total}] Failed to create triplet: {triplet_name}")
        except Exception as e:
            logger.error(f"[{idx + 1}/{total}] Error processing {triplet_name}: {e}")

    logger.info(f"Successfully generated {success_count}/{total} triplets")
    return success_count
