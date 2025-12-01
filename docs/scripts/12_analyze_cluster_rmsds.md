# 12_聚类内 RMSD 分析 (Cluster Intra-RMSD Analysis)

# 1. 目标与输入输出

## 目标
计算每个聚类（Cluster）内部样本之间的两两 RMSD（Root Mean Square Deviation），以量化评估聚类的结构一致性（紧密度）。通过统计数据，我们可以快速判断哪些簇包含高度相似的表位，哪些簇虽然被归为一类但结构差异较大（“相似的抗原有多不像”）。

## 输入
*   **聚类汇总表格**: `data/cluster/epitope_cluster_summary.csv` (由 `10_cluster_epitopes.py` 生成)
*   **表位结构数据 (NPZ)**: `data/extract_epitope/pockets/` (用于提取 Cα 坐标)。

## 输出
*   **总体统计表格**: `data/analysis/rmsd/cluster_rmsd_stats.csv`
    *   记录每个簇的宏观统计量：`cluster_id`, `size`, `mean_rmsd`, `std_rmsd`, `min_rmsd`, `max_rmsd`。
*   **簇内详细表格**: `data/analysis/rmsd/details/cluster_{id}.csv`
    *   记录簇内每一对样本的 RMSD：`epitope_1`, `epitope_2`, `rmsd`。

# 2. 技术方案

## 核心工具
*   **计算**: `numpy`, `biopython` (SVDSuperimposer)。
*   **并行化**: `joblib` (按簇并行计算)。

## 流程逻辑

### 1. 数据准备
*   读取 `epitope_cluster_summary.csv`。
*   筛选出 `size >= 2` 的簇。

### 2. 处理每个簇 (Per-Cluster Calculation)
对于每个簇：
1.  **加载坐标**: 读取簇内所有成员的 `.npz` 文件，提取 Cα 坐标。
2.  **检查兼容性**: 确保比较的表位具有相同数量的 Cα 原子（理论上聚类步骤已过滤，但需再次确认）。若数量不同，记录警告并跳过该对。
3.  **两两计算 (Pairwise Calculation)**:
    *   对于簇内所有唯一的配对 $(i, j)$ where $i < j$：
    *   计算 RMSD (使用 Kabsch 叠加)。
    *   收集结果：`[epitope_id_i, epitope_id_j, rmsd_value]`。
4.  **保存详细表格**: 将该簇的所有配对结果保存为 `details/cluster_{id}.csv`。
5.  **计算统计量**: 基于收集到的所有 RMSD 值，计算 mean, std, min, max。
6.  **返回统计结果**: 供主进程汇总。

### 3. 汇总与保存
*   收集所有簇的统计结果。
*   保存为 `cluster_rmsd_stats.csv`。

# 3. 进度追踪 (Implementation Progress)

- [ ] **环境准备**: 确认依赖。
- [ ] **主流程脚本 (`scripts/12_analyze_cluster_rmsds.py`)**:
    - [ ] 参数解析。
    - [ ] 实现数据加载与分组。
    - [ ] 实现 `process_cluster_rmsd` 函数。
    - [ ] 实现并行计算。
    - [ ] 结果汇总保存。
- [x] **测试与验证**:
    - [x] 运行脚本。
    - [x] 检查统计表格合理性。

# 4. 结果分析 (Results Analysis)

分析已完成，结果保存在 `data/analysis/rmsd/cluster_rmsd_stats.csv`。

**初步观察**:
*   **紧密簇 (Tight Clusters)**: 部分簇（如 Cluster 18）的 Mean RMSD 极低（~0.25 Å），说明其内部成员在表位结构上几乎完全一致，可能是同一抗原的高重复测定。
*   **多样簇 (Diverse Clusters)**: 部分簇（如 Cluster 9）的 RMSD 跨度较大（0.16 - 2.21 Å），说明虽然它们在聚类时满足了 RMSD < 2.5 Å 的阈值，但在细节上仍有显著构象差异，值得进一步通过对齐结构（Step 11 的产物）进行目视检查。
