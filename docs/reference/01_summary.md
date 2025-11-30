# summary.py 脚本说明

该脚本是数据预处理流程的核心部分，用于处理 SAbDab 的原始摘要文件（summary file）和对应的 PDB 结构文件。它负责数据清洗、序列提取、CDR 掩码标记以及格式标准化。

## 核心功能

1.  **数据过滤与加载**：
    *   读取 SAbDab 的 TSV 摘要文件。
    *   根据抗原类型（`ALLOWED_ANTIGEN_TYPES`）和分辨率（`RESOLUTION_THRESHOLD`，默认 4.5Å）过滤条目。
    *   剔除已知的错误 PDB ID（如 `4hjj`）。

2.  **结构解析 (Biopython)**：
    *   使用 `Bio.PDB` 解析 `.pdb` 文件。
    *   处理非标准氨基酸残基，将其映射为标准氨基酸（参考 `STANDARD_RESIDUE_SUBSTITUTIONS_INCASEOF_NON_STANDARD_RESIDUE`）。
    *   截断过长的序列以适配 Chothia 编号方案（VL <= 109, VH <= 113）。

3.  **序列处理与 CDR 掩码**：
    *   利用 `abnumber` 库对序列进行 Chothia 编号。
    *   生成掩码序列（masked sequence），将 CDR 区域（CDR1, CDR2, CDR3）替换为特定字符（默认为 'X'）。
    *   验证 CDR/FR 区域的完整性，不完整的条目将被丢弃。

4.  **错误修正**：
    *   **重链/轻链互换检测**：检测 SAbDab 数据中重链和轻链标注反转的情况，并自动修正（`exchange_heavy_and_light_info_for_single_entry`）。

5.  **并行处理与输出**：
    *   使用 `joblib` 并行处理每个条目。
    *   生成包含详细信息的 JSON 和 CSV 文件。
    *   记录处理失败的条目（如 PDB 解析失败、非标准残基异常等）。

## 主要函数

*   `parse_biopython_structure(...)`: 解析 PDB 实体链，提取氨基酸序列。
*   `mask_cdr(seq)`: 使用 `abnumber` 识别 CDR 区域并生成掩码序列。
*   `correct_sabdab(json_content)`: 检查并修正重轻链标注错误的条目。
*   `process_entry(...)`: 处理单个数据条目的主逻辑，包含结构加载、序列提取和初步验证。
*   `process_raw_data(args)`: 主流程控制，负责参数解析、并行任务分发和结果保存。

## 使用参数

*   `--chothia_dir`: 存放 Chothia 编号的 PDB 文件目录。
*   `--summary_dir`: SAbDab 摘要文件路径（TSV格式）。
*   `--output_dir`: 输出目录。
*   `--output_file_name`: 输出文件名（默认为 `summary.json`）。
*   `--consider_no_antigen`: 标志位，如果设置则包含没有抗原的条目。

## 输出

脚本执行后，将在指定的输出目录（`--output_dir`）生成以下文件结构：

### 文件结构

```text
output_dir/
├── summary.json                 # 处理后的核心数据文件（JSON格式）
├── summary.csv                  # 对应的表格文件（CSV格式），方便查看
├── failed_to_process_with_bio.json # 处理失败的条目名称列表
└── logs-for-preprocess/         # 日志文件夹
    └── processing_YYYYMMDD_HHMMSS.log # 详细的运行日志
```

### JSON/CSV 字段说明

`summary.json` 中的每个条目（Key 为 `entry_name`）和 `summary.csv` 的每一行包含以下字段：

*   **`pdb`**: PDB ID（小写）。
*   **`file_name`**: (仅 CSV) 条目的唯一标识符，格式为 `{Pdb_id}_{H_chain}_{L_chain}_{Antigen_chains}`。
*   **`H_chain_id`**: 重链的链 ID。
*   **`L_chain_id`**: 轻链的链 ID。
*   **`H_chain_seq`**: 提取出的重链氨基酸序列（基于 Chothia 编号，去除了非 Fv 区域）。
*   **`L_chain_seq`**: 提取出的轻链氨基酸序列（基于 Chothia 编号，去除了非 Fv 区域）。
*   **`H_chain_masked_seq`**: 掩码处理后的重链序列（CDR 区域被替换为 'X'）。
*   **`L_chain_masked_seq`**: 掩码处理后的轻链序列（CDR 区域被替换为 'X'）。
*   **`antigen_chain_id`**: 抗原链 ID 的列表（List）。
*   **`antigen_seq`**: (JSON 中) 抗原链序列的字典，Key 为链 ID，Value 为序列；(CSV 中) 序列的字符串表示。
*   **`antigen_type`**: 抗原类型（例如 "protein"）。
*   **`resolution`**: 晶体结构的分辨率（单位：Å）。
*   **`scfv`**: 布尔值，指示是否为 scFv 结构。
*   **`date`**: 结构发布的日期。
*   **`index_in_summary`**: 该条目在原始 SAbDab 摘要文件中的索引位置。

### 错误处理

*   **`failed_to_process_with_bio.json`**: 记录了那些因为 PDB 文件解析失败、文件缺失或其他 Biopython 处理错误而被丢弃的条目。
*   **日志文件**: 记录了每个处理失败条目的具体原因（例如 "Heavy chain error... No module named 'anarci'"）。
