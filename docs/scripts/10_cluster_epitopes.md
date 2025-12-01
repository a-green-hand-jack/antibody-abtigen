# 10_抗原表位聚类与对齐 (Epitope Clustering & Alignment)

# 1. 目标与输入输出

## 目标
对已提取的抗原表位进行聚类，识别结构和序列相似的表位组，以便后续的代表性选择、去冗余和数据集构建。

## 输入
*   **表位汇总表格**: `data/extract_epitope/epitope_summary.csv` (由 `09_extract_epitope.py` 生成)
    *   包含 `epitope_id`, `pocket_npz_path` 等关键信息。
*   **结构数据 (NPZ)**: `data/extract_epitope/pockets/{epitope_id}.npz`
    *   存储裁剪后的原子级结构信息 (coords, atom_types, etc.)。

## 输出
*   **聚类后的汇总表格**: `data/cluster/epitope_cluster_summary.csv`
    *   在原始汇总表格的基础上增加 `cluster_id`, `is_representative` 等字段。
*   **(可选) 聚类统计报告**: 记录每个簇的大小、代表性选择标准等。

# 2. 技术方案

## 核心工具
*   **数据处理**: `pandas`, `numpy`
*   **结构比对**: `biopython` (Kabsch 算法 for RMSD), 或自定义实现。
*   **序列比对**: `Bio.pairwise2` 或 `skbio.sequence.distance` (用于计算序列同一性)。
*   **图算法**: `networkx` (用于构建图和查找连通分量)。
*   **并行计算**: `joblib` (用于加速相似度矩阵计算)。

## 流程逻辑

### 1. 读取数据
*   加载 `data/extract_epitope/epitope_summary.csv`。
*   对于每个表位，读取其对应的 `.npz` 文件以获取结构数据。

### 2. 构建相似度矩阵 (Pairwise Similarity Calculation)
*   **粗筛 (Residue Count Pre-filtering)**:
    *   为减少计算量，首先根据表位中残基数量的差异进行初步过滤。
    *   只有当两个表位的残基数量差值在设定阈值内 (例如，`abs(res_count1 - res_count2) <= N_res_threshold`) 时，才进行详细的结构和序列比较。
    *   **注意**: 这里的残基数量指 `09_extract_epitope.py` 步骤中提取并过滤后的表位残基数量。
*   **结构比对 (Structural Alignment)**:
    *   对于通过粗筛的表位对，提取它们的 Cα 原子坐标。
    *   使用 Kabsch 算法计算最佳刚体对齐后的 **RMSD (Root Mean Square Deviation)**。
    *   **阈值**: 仅考虑 `RMSD <= 2.5 Å` (可配置) 的表位对。
*   **序列比对 (Sequence Alignment)**:
    *   对于通过粗筛的表位对，提取其氨基酸序列。
    *   计算两个表位序列的 **Sequence Identity**。
    *   **阈值**: 仅考虑 `Sequence Identity >= 40%` (可配置) 的表位对。
*   **并行化**: 使用 `joblib` 或类似库并行计算相似度矩阵。

### 3. 构建相似性图 (Similarity Graph Construction)
*   **节点**: 每个独特的表位 (`epitope_id`) 作为一个图的节点。
*   **边**: 如果两个表位同时满足以下条件，则在它们之间添加一条边：
    *   `RMSD <= RMSD_threshold`
    *   `Sequence Identity >= SeqID_threshold`

### 4. 聚类 (Clustering)
*   使用 `networkx` 库或其他图算法，通过查找图的 **连通分量 (Connected Components)** 来进行聚类。每个连通分量代表一个表位簇。
*   为每个簇分配一个唯一的 `cluster_id`。

### 5. 选择簇代表 (Representative Selection)
*   每个簇中选择一个代表性表位。选择标准可能包括：
    *   簇内残基数量最接近中位数/平均值的表位。
    *   原始数据集中质量评分较高（如果有的话）的表位。
    *   简单起见，可以先选择 `pdb_id` 字典序最小的表位作为代表。

### 6. 输出结果
*   将 `cluster_id` 和 `is_representative` 信息添加回 `epitope_summary.csv`，并保存为 `data/cluster/epitope_cluster_summary.csv`。

## 阈值设定 (Proposed Defaults, configurable)
*   `N_res_threshold` (for coarse pre-filtering): 5 residues
*   `RMSD_threshold`: 2.5 Å
*   `SeqID_threshold`: 40%

# 3. 进度追踪 (Implementation Progress)

- [x] **环境准备**: 确认 `biopython`, `networkx` 依赖。
- [x] **主流程脚本 (`scripts/10_cluster_epitopes.py`)**:
    - [x] 参数解析 (输入输出路径, 阈值配置)。
    - [x] 实现数据加载。
    - [x] 实现表位残基数量粗筛逻辑。
    - [x] 实现结构比对 (Kabsch, RMSD)。
    - [x] 实现序列比对 (Sequence Identity)。
    - [x] 实现相似度矩阵并行计算。
    - [x] 实现相似性图构建。
    - [x] 实现连通分量聚类。
    - [x] 实现簇代表选择。
    - [x] 结果汇总与 CSV 输出。
- [x] **文档 (`docs/scripts/10_cluster_epitopes.md`)**:
    - [x] 初稿完成并与用户确认。
- [x] **测试与验证**:
    - [x] 运行小规模测试集合。
    - [x] 验证聚类结果的合理性。
    - [x] 检查日志。

# 4. 结果分析 (Results Analysis)

截至 2025年12月1日，我们对 **3,537** 个抗原表位进行了聚类分析。
使用的阈值参数为：RMSD ≤ 2.5 Å, Sequence Identity ≥ 40%, Residue Count Diff ≤ 5。

**聚类概况**:
*   **总聚类数**: 2,895 个
*   **单体聚类 (Singletons)**: 2,558 个 (约占总数据的 72%，显示出数据集中具有较高的结构/序列多样性)。
*   **相似表位组 (Clusters with >1 item)**: **337** 个。这些组共包含 **979** 个抗原表位。

**聚类大小分布 (Cluster Size Distribution)**:
| 每个 Group 内的表位数量 (Size) | 这样的 Group 有多少个 (Count) | 总涉及表位数量 | 说明 |
| :--- | :--- | :--- | :--- |
| **2** | **235** | 470 | 成对出现的相似表位 |
| **3** | **52** | 156 | |
| **4** | **24** | 96 | |
| **5** | **9** | 45 | |
| **6** | **4** | 24 | |
| **7** | **3** | 21 | |
| **8** | **1** | 8 | |
| **9** | **3** | 27 | |
| **10** | **1** | 10 | |
| **13** | **1** | 13 | |
| **18** | **1** | 18 | |
| **26** | **1** | 26 | 潜在的热门抗原或重复测定 |
| **27** | **1** | 27 | |
| **38** | **1** | 38 | 最大簇 |