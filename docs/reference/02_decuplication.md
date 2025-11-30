# decuplication.py 脚本说明

该脚本主要用于对 SAbDab 数据集中的抗体条目进行**去重**处理。它基于抗体的互补决定区（CDR）序列的相似性，移除那些在同一 PDB ID 下序列过于相似的冗余条目。

## 核心功能

1.  **基于 CDR 的相似度计算**：
    *   提取抗体重链的 CDR1, CDR2, CDR3 序列。
    *   将三个 CDR 拼接成一个长序列。
    *   使用 **Levenshtein 距离**（编辑距离）计算两个序列之间的差异。

2.  **去重逻辑**：
    *   按 `H_chain_seq`（重链全序列）长度降序排列同一 PDB ID 下的所有条目。
    *   保留最长的序列作为初始参考。
    *   遍历剩余条目，计算其拼接 CDR 序列与已保留条目的 CDR 序列之间的 Levenshtein 距离。
    *   如果距离小于或等于设定的阈值（`distance_threshold`），则视为重复并丢弃；否则保留该条目。

3.  **并行处理**：
    *   脚本使用 `joblib` 对不同的 PDB ID 进行并行处理，提高大规模数据处理的效率。

4.  **日志记录**：
    *   详细记录了被保留的条目以及因相似度过高而被丢弃的条目对比信息（包括序列比对结果）。

## 主要函数

*   `levenshtein_distance(seq1, seq2)`: 计算两个字符串之间的编辑距离。
*   `get_cdr(original_seq, masked_seq)`: 根据原始序列和掩码序列（CDR 位置标记为 'X'）提取 CDR1, CDR2, CDR3。
*   `decuplicate(pdb_id, group_df, distance_threshold)`: 处理单个 PDB ID 下的数据组，执行具体的去重逻辑。
*   `main(args)`: 主程序入口，负责加载数据、并行调度、汇总结果并保存到 CSV 文件。

## 使用参数

*   `--distance_threshold`: 判定重复的 Levenshtein 距离阈值（默认为 9）。
*   `--original_data_path`: 原始 CSV 数据文件路径（默认为 `./data/summary.csv`）。
*   `--output_dir`: 输出目录。
*   `--output_file_name`: 输出文件名（不带后缀）。

## 输出

*   去重后的 CSV 文件（例如 `summary-decuplication-distance_threshold_9.csv`）。
*   详细的日志文件，位于 `logs-for-decuplication/` 目录下。
