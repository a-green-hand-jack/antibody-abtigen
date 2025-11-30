# 04_cutoff.py 脚本说明

该脚本用于根据时间截止日期（Cutoff Date）筛选 PDB 条目，这对于防止机器学习模型训练中的**数据泄露**（Data Leakage）至关重要。它通过检查 PDB 文件的沉积（deposition）和发布（release）日期来实现筛选。
这个脚本的作用不是过滤而是为了分割数据集为训练、验证、测试数据集做准备。

## 核心功能

1.  **日期提取**：
    *   使用 `gemmi` 库解析 `.cif` 格式的结构文件。
    *   从 CIF 文件的 `_pdbx_database_status` 和 `_pdbx_audit_revision_history` 字段中提取沉积日期、发布日期和最后修订日期。

2.  **时间筛选**：
    *   读取 SAbDab 摘要文件中的日期信息。
    *   将摘要中的日期以及 CIF 文件中的发布日期与用户指定的 `cutoff_date` 进行比较。
    *   分别记录在 SAbDab 记录中早于截止日期的条目，以及在结构文件元数据中早于截止日期的条目。

3.  **一致性检查**：
    *   比较基于 SAbDab 日期和基于 CIF 内部日期的筛选结果是否一致（即 Boltz 筛选集是否是 SAbDab 筛选集的子集）。

4.  **结果保存**：
    *   将发布日期在截止日期**之前（包含当天）**的 PDB ID 列表保存为 JSON 文件 (`before_cutoff_in_sabdab.json`)。这意味着该文件记录的是可用于训练模型的数据（即 cutoff date 之前的数据）。

## 主要函数

*   `get_dates(block)`: 从 gemmi 解析的 CIF block 中提取日期信息。
*   `get_pdb_id_before_cutoff_date(...)`: 遍历所有条目，对比日期并返回符合条件的 ID 集合。
*   `main(args)`: 主程序入口。

## 使用参数

*   `--summary_file_path`: 摘要 CSV 文件路径。
*   `--cif_dir`: **CIF 文件** 所在的目录（注意：此脚本依赖 `.cif` 文件而非 `.pdb`）。
*   `--out_dir`: 输出目录。
*   `--cutoff_date`: 截止日期字符串（格式 `YYYY-MM-DD`，默认 `2021-09-30`）。

## 依赖
*   `gemmi`: 用于高效解析 CIF 文件。
