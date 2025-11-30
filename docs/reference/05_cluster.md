# 05_cluster.py 脚本说明

该脚本利用 **MMseqs2** 工具对抗体序列（特别是 CDR3 区域）进行聚类。聚类是构建非冗余数据集的关键步骤，有助于减少模型训练偏差。

## 核心功能

1.  **提取 CDR 序列**：
    *   从处理后的 CSV 文件中读取重链和轻链序列。
    *   根据掩码序列（masked sequence）解析出 CDR1, CDR2, CDR3 片段。

2.  **生成 FASTA 文件**：
    *   将提取出的各 CDR 序列分别写入独立的 FASTA 文件（如 `heavy_chain_cdr3.fasta`）。
    *   FASTA 头部包含文件名和长度信息。

3.  **调用 MMseqs2 聚类**：
    *   使用 Python `subprocess` 调用外部命令 `mmseqs easy-cluster`。
    *   根据设定的序列一致性阈值（`seq_identity_threshold`）对 FASTA 文件进行聚类。

## 主要函数

*   `get_cdr(original_seq, masked_seq)`: 从全序列中切分出 CDR 区域。
*   `processed_csv_to_fasta_only_cdr3(...)`: 生成包含 CDR 序列的 FASTA 文件。
*   `use_mmseqs_to_cluster(...)`: 包装 MMseqs2 命令执行聚类操作。

## 使用参数

*   `--processed_csv_path`: 输入的 CSV 数据文件路径。
*   `--decuplication_threshold`: 去重阈值（用于校验文件名一致性）。
*   `--out_dir`: 输出目录。
*   `--seq_identity_threshold`: 序列一致性阈值（默认 `0.5`）。
*   `--mmseqs_path`: MMseqs2 可执行文件的路径（默认 `./mmseqs/bin/mmseqs`）。

## 输出
*   在输出目录下生成按去重阈值命名的文件夹。
*   包含 CDR 序列的 FASTA 文件。
*   MMseqs2 的聚类结果（`.tsv` 格式），包含聚类中心和成员信息。
