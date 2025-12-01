# 13_结构清洗与拆分 (Structure Cleaning & Splitting)

# 1. 目标与输入输出

## 目标
对已对齐的抗原-抗体复合物结构进行清洗和拆分。
1.  **清洗**: 去除水分子、溶剂和无关组件。
2.  **筛选**: 剔除不参与相互作用的“非抗原链”和无关蛋白。
3.  **拆分**: 将抗体和抗原分离为独立的文件，但**严格保持**其在对齐步骤中确定的相对空间位置。

## 输入
*   **对齐后的结构文件**: `data/aligned/{cluster_id}/{epitope_id}.cif` (由 `11_align_structures.py` 生成)。
*   **元数据**: `data/cluster/epitope_cluster_summary.csv` (提供 H/L 链和抗原链的 ID)。

## 输出
*   **目录结构**: `data/cleaned_split/{cluster_id}/{epitope_id}/`
*   **文件**:
    *   `antibody.cif`: 仅包含抗体链 (H, L)。
    *   `antigen.cif`: 仅包含指定的抗原链。

# 2. 技术方案

## 核心工具
*   **Gemmi**: 用于结构读取、链筛选、去除水分子和文件保存。

## 流程逻辑

### 1. 数据加载
*   读取 `epitope_cluster_summary.csv`。
*   遍历每一行数据。

### 2. 处理流程 (Per-Epitope Processing)
对于每个表位条目：
1.  **路径解析**:
    *   输入路径: `data/aligned/{cluster_id}/{epitope_id}.cif`。
    *   如果输入文件不存在（可能该簇太小被跳过），则跳过。
2.  **解析链信息**:
    *   获取 `H_chain`, `L_chain`。
    *   解析 `antigen_chains` (处理 "A|B" 格式)。
3.  **结构读取与清洗**:
    *   使用 `gemmi.read_structure` 读取 CIF。
    *   调用 `structure.remove_waters()` 移除水分子。
    *   (可选) 移除氢原子以减小体积。
4.  **拆分构建**:
    *   创建两个新的 `gemmi.Structure` 对象: `st_ab` (Antibody) 和 `st_ag` (Antigen)。
    *   遍历原结构中的所有 Chain：
        *   如果 Chain ID 匹配 H 或 L -> `clone` 并加入 `st_ab`。
        *   如果 Chain ID 在抗原链列表中 -> `clone` 并加入 `st_ag`。
        *   其他链 -> 忽略。
5.  **保存**:
    *   创建输出文件夹 `data/cleaned_split/{cluster_id}/{epitope_id}/`。
    *   保存 `st_ab` 为 `antibody.cif`。
    *   保存 `st_ag` 为 `antigen.cif`。

# 3. 进度追踪 (Implementation Progress)

- [x] **环境准备**: 确认 `gemmi`, `pandas` 依赖。
- [x] **主流程脚本 (`scripts/13_clean_and_split.py`)**:
    - [x] 参数解析。
    - [x] 实现链 ID 解析逻辑。
    - [x] 实现 Gemmi 结构拆分逻辑 (Clone Chains)。
    - [x] 实现并行处理。
- [x] **测试与验证**:
    - [x] 运行脚本。
    - [x] 检查输出目录结构。

# 4. 问题排查 (Troubleshooting)

*   **空文件问题 (Empty CIFs)**:
    *   **现象**: 生成的 CIF 文件仅几百字节，只包含 header 无原子。
    *   **原因**: `gemmi.Structure.add_model(model)` 会创建模型的**副本**。如果随后向原始 `model` 变量添加链，这些链不会进入 Structure 中。
    *   **解决**: 在添加 Model 后，通过 `model = structure[0]` 获取 Structure 内部的引用，再进行操作。此外，写入前必须调用 `setup_entities()` 以生成合法的 mmCIF。

# 5. 验证 (Verification)

脚本运行完成后，检查 `data/cleaned_split` 目录。以 Cluster 34 为例，应包含 38 个子文件夹，每个文件夹内包含 `antigen.cif` 和 `antibody.cif`。

```bash
ls -R data/cleaned_split/34 | grep ".cif" | wc -l
# Output should be 76 (38 epitopes * 2 files)
```

这些文件已经去除了水分子，并且只包含特定的抗原/抗体链。它们的坐标保持了 Step 11 中的对齐状态。
