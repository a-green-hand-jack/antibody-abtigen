# 11_结构对齐与保存 (Structure Alignment & Saving)

# 1. 目标与输入输出

## 目标
基于聚类结果，将同一簇（Cluster）内的所有抗原-抗体复合物根据其**抗原表位**的结构进行刚体对齐，并将对齐后的完整结构保存为 CIF 文件，以便于可视化分析（如 PyMOL）和后续的结构比较。

## 输入
*   **聚类汇总表格**: `data/cluster/epitope_cluster_summary.csv` (由 `10_cluster_epitopes.py` 生成)
    *   关键字段: `cluster_id`, `is_representative`, `pdb_path`, `pocket_npz_path`, `epitope_id`.
*   **原始结构文件**: CIF 格式 (路径由 `pdb_path` 指定)。
*   **表位结构数据**: NPZ 格式 (路径由 `pocket_npz_path` 指定，用于提取 Cα 坐标计算对齐矩阵)。

## 输出
*   **对齐后的结构文件**: 
    *   目录结构: `data/aligned/{cluster_id}/{epitope_id}.cif`
    *   内容: 经过旋转平移的完整抗原-抗体复合物结构。

# 2. 技术方案

## 核心工具
*   **结构处理**: `gemmi` (用于读取 CIF, 应用变换矩阵, 保存 CIF)。
*   **对齐算法**: `biopython` (`Bio.SVDSuperimposer`) 或 `numpy` (实现 Kabsch 算法)。
*   **并行计算**: `joblib` (按簇并行处理)。

## 流程逻辑

### 1. 数据准备
*   读取 `epitope_cluster_summary.csv`。
*   过滤掉 `cluster_id == -1` (未聚类) 的条目（如果有）。
*   按 `cluster_id` 进行分组。

### 2. 处理每个簇 (Per-Cluster Processing)
对于每个簇：
1.  **创建输出目录**: `data/aligned/{cluster_id}/`.
2.  **识别参考结构 (Reference)**:
    *   找到 `is_representative == True` 的条目。
    *   读取其 `.npz` 文件获取表位 Cα 坐标 (`Ref_Coords`)。
    *   读取其原始 `.cif` 文件。
    *   将原始 CIF 直接保存到输出目录（无需变换）。
3.  **对齐其他成员 (Mobile)**:
    *   对于簇内的每个非代表成员：
        *   读取其 `.npz` 文件获取表位 Cα 坐标 (`Mobile_Coords`)。
        *   **检查原子一致性**: 确保 `Mobile_Coords` 和 `Ref_Coords` 的原子数量一致（基于之前的聚类逻辑，这应该已基本满足，但需防范异常）。
        *   **计算变换**: 使用 SVD/Kabsch 算法计算 `Mobile_Coords` 到 `Ref_Coords` 的旋转矩阵 ($R$) 和平移向量 ($t$)。
        *   **应用变换**: 读取该成员的原始 `.cif` 文件，对其中**所有原子**的坐标 $(x, y, z)$ 应用变换：$v' = R \cdot v + t$。
        *   **保存**: 将变换后的结构保存为 `data/aligned/{cluster_id}/{epitope_id}.cif`。

### 3. 优化策略
*   **最小簇大小**: 默认只处理 `Size >= 2` 的簇（单体簇不需要对齐，除非为了统一归档）。可通过参数 `--min_size` 控制。
*   **并行化**: 以“簇”为单位进行并行处理。

# 3. 进度追踪 (Implementation Progress)

- [x] **环境准备**: 确认 `gemmi`, `biopython` 依赖。
- [x] **主流程脚本 (`scripts/11_align_structures.py`)**:
    - [x] 参数解析 (输入输出路径, 最小簇大小)。
    - [x] 实现数据加载与分组。
    - [x] 实现对齐核心函数 (`align_and_save_cluster`)。
        - [x] 读取 NPZ 获取 Cα。
        - [x] 计算 RMSD/Transformation (SVD)。
        - [x] 使用 Gemmi 应用变换。
        - [x] 保存 CIF。
    - [x] 实现并行处理。
- [x] **测试与验证**:
    - [x] 运行脚本。
    - [x] 验证输出文件数量与聚类大小一致。

# 4. 验证 (Verification)

脚本运行完成后，可以检查 `data/aligned` 目录。例如，对于最大的簇（ID 34），应包含 38 个 CIF 文件：

```bash
ls data/aligned/34 | wc -l
# Output should be 38
```

你可以使用 PyMOL 打开该目录下的所有文件来直观验证对齐效果：
```bash
pymol data/aligned/34/*.cif
```
在 PyMOL 中，你应该能看到所有抗原的表位部分重合在一起，而抗体则根据其结合角度呈现出不同的空间排布。
