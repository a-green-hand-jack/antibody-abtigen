# download_sabdab.py 脚本说明

该脚本用于自动下载 SAbDab (Structural Antibody Database) 的摘要文件以及相关的抗体-抗原复合物 PDB 结构文件。它包含自动筛选逻辑，仅下载具有蛋白类抗原的结构。

## 核心功能

1.  **摘要文件下载**：
    *   检查本地是否存在 `data/meta/sabdab_summary_all.tsv`。
    *   如果不存在，从 OPIG 服务器下载最新的摘要文件。

2.  **数据筛选**：
    *   读取摘要文件。
    *   **筛选条件**：
        *   `antigen_type` 必须为 `"protein"` (仅保留蛋白抗原)。
        *   `antigen_chain` 必须存在且非空 (排除无抗原结合的单抗体结构)。
    *   提取符合条件的唯一 PDB ID 列表。

3.  **结构文件批量下载**：
    *   遍历筛选出的 PDB ID。
    *   从 SAbDab 服务器构建下载链接：`https://opig.stats.ox.ac.uk/webapps/abdb/entries/{pdb_id}/structure/{pdb_id}.pdb`。
        *   *注意*：脚本下载的是 **Chothia 编号** 格式的 PDB 文件。
    *   将文件保存到 `data/raw_data/chothia/` 目录下（文件名格式为 `{pdb_id}.pdb`）。
    *   支持断点续传（跳过已存在的文件）。

## 主要函数

*   `download_file(url, filename)`: 通用下载辅助函数，包含进度条显示（`tqdm`）、错误处理和不完整文件清理功能。
*   `main()`: 主流程控制，执行摘要下载、数据清洗筛选和 PDB 批量下载。

## 输出目录结构

脚本运行时会自动创建以下目录（如果不存在）：

*   `data/meta/`: 存放 `sabdab_summary_all.tsv`。
*   `data/raw_data/chothia/`: 存放下载的 Chothia 编号 `.pdb` 文件。

## 使用方法

直接运行脚本即可：

```bash
python scripts/download_sabdab.py
```

## 依赖库

*   `requests`: 处理 HTTP 请求。
*   `pandas`: 读取 TSV 和筛选数据。
*   `tqdm`: 显示下载进度条。
