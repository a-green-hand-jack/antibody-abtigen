# Python 项目执行计划：人源-鼠源抗原结构配对数据集构建

根据要求，我们将数据处理流程划分为五个阶段。每个阶段说明**输入/输出**、**处理逻辑**、**使用的 Python 模块与第三方库**、**文件命名及目录结构**，并阐述如何实现**增量构建机制**。

> **约束条件**：
> - 全程在服务器环境下运行
> - 仅使用公开数据库和已有模型（如 SAbDab、AlphaFold DB）
> - 不调用任何额外的 AI 建模工具
> - 所有结构数据统一使用 **mmCIF 格式**处理和存储（**不使用 PDB 格式**）

---

## 阶段一：数据下载与原始数据存储

### 输入及数据源
- **无本地输入**，通过网络从以下数据库获取：
  - **SAbDab**（Structural Antibody Database）：
    - 提供抗原-抗体复合物的元数据（TSV/CSV）
    - 从 SAbDab 获取复合物清单，包含 PDB 编号、链信息、抗原名称、物种等
    - 通过 RCSB PDB 或 PDBe 的公开接口下载对应 mmCIF 坐标文件（如 `https://files.rcsb.org/download/{PDB_ID}.cif`）

### 处理逻辑
1. **获取清单**：使用 `requests` 下载 SAbDab 提供的汇总 TSV 文件。
2. **批量下载 mmCIF**：
   - 用 `pandas` 读取 TSV，提取 PDB 编号列表
   - 对每个 PDB 编号，构造 URL 下载 `.cif` 文件并保存为 `{PDB_ID}.cif`
3. **保存元数据**：
   - 保存原始 TSV 为 `sabdab_summary.tsv`
   - 可选：提取关键字段（如 PDB ID、抗原链、抗原物种、抗体链）生成 `sabdab_summary_filtered.csv`

### 输出
```
./data/antigen_antibody/SAbDab/raw/
├── sabdab_summary.tsv
├── sabdab_summary_filtered.csv
├── 1XYZ.cif
├── 1ABC.cif
└── ...
```

### 文件命名规范
- 结构文件：`{PDB_ID}.cif`（如 `1XYZ.cif`）
- 元数据：`sabdab_summary.tsv`、`sabdab_summary_filtered.csv`

### 使用模块与库
- `requests`：HTTP 下载
- `pandas`：元数据处理
- `os` / `pathlib`：路径与目录管理

### 增量构建机制
- 下载每个 `.cif` 前，检查是否已存在，若存在则跳过
- 可记录上次下载时间或文件哈希，仅更新新增条目
- 维护“已处理 PDB 列表”，仅下载新增复合物

---

## 阶段二：数据清理

### 输入
- `./SAbDab/raw/` 下的 mmCIF 文件
- 对应的元数据表（TSV/CSV），提供抗原/抗体链 ID

### 处理逻辑
1. **解析结构**：使用 `Bio.PDB.MMCIFParser` 读取 mmCIF
2. **确定保留链**：从元数据中提取抗原链、抗体（重链/轻链）链 ID
3. **过滤链**：
   - 仅保留目标链（抗原 + 抗体）
   - 移除水分子（`HOH`）、配体、非蛋白残基（`HETATM`）
4. **保存清理结构**：用 `Bio.PDB.MMCIFIO` 写出精简后 mmCIF

### 输出
```
./data/antigen_antibody/SAbDab/cleaned/
├── 1XYZ.cif   # 仅含抗原+抗体链
├── 1ABC.cif
└── ...
```

### 文件命名规范
- 保留原始 PDB 编号：`1XYZ.cif`

### 使用模块与库
- `Bio.PDB`（`MMCIFParser`, `MMCIFIO`）
- 或 `gemmi`（支持 mmCIF 的高效结构处理）
- `os` / `pathlib`

### 增量构建机制
- 若 `cleaned/1XYZ.cif` 已存在，且清理逻辑未变，则跳过
- 可通过时间戳或文件对比判断是否需重处理
- 扫描 `raw/` 与 `cleaned/` 目录差异，仅处理缺失项

---

## 阶段三：拆分抗原与抗体

### 输入
- `./SAbDab/cleaned/` 下的复合物 mmCIF 文件

### 处理逻辑
1. **加载结构**：解析每个清理后的 mmCIF
2. **提取抗原链**：根据元数据获取抗原链 ID（可能为多条链）
3. **提取抗体链**：通常包含一条重链 + 一条轻链
4. **分别保存**：
   - 抗原 → `antigen/{PDB_ID}_antigen.cif`
   - 抗体 → `antibody/{PDB_ID}_antibody.cif`

### 输出
```
./data/antigen_antibody/SAbDab/antigen/
├── 1XYZ_antigen.cif
├── 1ABC_antigen.cif
└── ...

./data/antigen_antibody/SAbDab/antibody/
├── 1XYZ_antibody.cif
├── 1ABC_antibody.cif
└── ...
```

### 文件命名规范
- `{PDB_ID}_antigen.cif`
- `{PDB_ID}_antibody.cif`

### 使用模块与库
- `Bio.PDB`（支持链级提取与保存）
- 可选：`MDAnalysis` 或 `PyMOL`（需支持 mmCIF 输出）

### 增量构建机制
- 检查 `antigen/` 和 `antibody/` 是否已存在对应文件
- 仅缺失部分重新生成
- 例如：若 `1XYZ_antigen.cif` 缺失，但抗体存在，仅提取抗原

---

## 阶段四：人源-鼠源抗原对齐

### 输入
- **人源抗原结构**：来自 `antigen/` 且元数据中标注为 **Human**
- **鼠源抗原结构**：通过 **AlphaFold DB**（或 PDB）获取

### 处理逻辑
1. **识别鼠源同源蛋白**：
   - 若人源抗原有 UniProt ID，调用 UniProt ID mapping 获取小鼠同源 ID
   - 否则基于抗原名称（如 "IL2" → "MouseIL2"）映射
2. **下载鼠源结构**：
   - 从 AlphaFold DB 下载 mmCIF：
     `https://alphafold.ebi.ac.uk/files/AF-{UniProtID}-F1-model_v4.cif`
   - 可使用 `Bio.PDB.alphafold_db` 或 `requests`
3. **结构对齐**：
   - 提取人源 & 鼠源序列（`get_sequence()`）
   - 序列比对（`Bio.Align`）
   - 基于 Cα 坐标使用 `Superimposer` 计算旋转/平移
   - 应用变换到鼠源结构
   - （可选）使用 `PyMOL` 的 `cmd.align()` 简化流程
4. **保存配对结果**：
   - 为每对创建子目录：`{human_ag}_{chain}_{mouse_ag}`
   - 存放两个文件：
     - `human_ag_chain.cif`
     - `mouse_ag_chain_aligned.cif`
   - 更新汇总 CSV：`human_mouse_pairs.csv`

### 输出
```
./data/antigen_antibody/HumanMouse/
└── IL2_A_MouseIL2/
    ├── human_ag_chain.cif
    └── mouse_ag_chain_aligned.cif
```

并生成：
- `human_mouse_pairs.csv`（含人源PDB、链、鼠源名称、来源等）

### 使用模块与库
- `Bio.PDB`（结构解析、序列提取、Superimposer）
- `Bio.Align`（序列比对）
- `requests`（下载 AlphaFold 模型）
- `pandas`（CSV 记录）
- `os` / `pathlib`
- 可选：`pymol2`（PyMOL 无头模式对齐）

### 增量构建机制
- 若 `{human}_{chain}_{mouse}/` 目录已存在，跳过
- 鼠源模型下载后缓存，避免重复请求
- 新增人源抗原会触发新对齐任务
- CSV 表追加新记录，不覆盖已有

---

## 阶段五：构建完整的鼠源抗原结构

### 输入
- `HumanMouse/` 下的所有对齐结果
- 特别关注多链抗原（如由链 A + B 组成）

### 处理逻辑
1. **按抗原名称分组**：
   - 例如 `IL2_A_MouseIL2` 和 `IL2_B_MouseIL2` 属于同一抗原
2. **拼接鼠源链**：
   - 读取各 `mouse_ag_chain_aligned.cif`
   - 合并为单一 `Structure` 对象（保留原始链 ID）
   - 若为独立亚基，直接合并；若为连续肽链，检查末端距离（可选）
3. **保存完整结构**：
   - 单链抗原：直接复制对齐文件为完整结构
   - 多链抗原：合并后保存为 `{mouse_ag_name}.cif`

### 输出
```
./data/antigen_antibody/MouseAntigen/
├── MouseIL2.cif
├── MouseAntigenX.cif
└── ...
```

### 使用模块与库
- `Bio.PDB`（结构合并：`structure[0].add(chain)`）
- 或 `gemmi`（支持多链模型）
- `itertools.groupby`（按抗原名分组）
- 日志记录：成功/失败抗原列表

### 增量构建机制
- 若 `MouseAntigen/{name}.cif` 已存在，跳过
- 若 `HumanMouse/` 新增链（如因新 PDB 加入），则重新构建该抗原
- 可通过比较链数量或维护链计数表判断是否需更新

---

## 文件目录结构总览

```
./data/antigen_antibody/
├── SAbDab/
│   ├── raw/
│   │   ├── sabdab_summary.tsv
│   │   ├── sabdab_summary_filtered.csv
│   │   ├── 1XYZ.cif
│   │   └── ...
│   ├── cleaned/
│   │   ├── 1XYZ.cif
│   │   └── ...
│   ├── antigen/
│   │   ├── 1XYZ_antigen.cif
│   │   └── ...
│   └── antibody/
│       ├── 1XYZ_antibody.cif
│       └── ...
├── HumanMouse/
│   └── IL2_A_MouseIL2/
│       ├── human_ag_chain.cif
│       └── mouse_ag_chain_aligned.cif
└── MouseAntigen/
    ├── MouseIL2.cif
    └── ...
```

**说明文件**：
- `sabdab_summary.tsv`（原始元数据）
- `sabdab_summary_filtered.csv`（筛选后）
- `human_mouse_pairs.csv`（人-鼠配对清单）

---

## 增量构建机制汇总

整个流程支持**幂等性**与**增量更新**：

| 阶段 | 增量策略 |
|------|--------|
| 阶段一 | 检查文件是否存在；仅下载新增 PDB |
| 阶段二 | 仅处理 `raw/` 中未清理的 PDB |
| 阶段三 | 仅生成缺失的 `_antigen.cif` 或 `_antibody.cif` |
| 阶段四 | 跳过已存在的配对目录；复用已下载鼠源模型 |
| 阶段五 | 按抗原名检查完整结构是否存在；支持链更新后重建 |

> ✅ **优势**：可定期运行完整管道，自动处理新增数据，避免重复计算，保证效率与一致性。




# 拓展模块计划

根据要求，我们将把数据处理流程划分为五个阶段，每个阶段说明输入/输出、处理逻辑、使用的Python模块和第三方库、文件命名及目录结构，并阐述如何实现增量构建机制。整个流程均在服务器环境下运行，仅使用公开数据库和已有模型（如 SAbDab、AlphaFold DB），不调用任何额外的AI建模工具。所有结构数据统一使用 **mmCIF** 格式处理和存储，不使用PDB格式。

此外，新增两个关键模块：

1. **人-鼠抗原表位一致性验证模块**（用于评估结构配对是否真实对应同一表位）；
2. **SAbDab 鼠源抗原复合物识别模块**（优先利用已存在的鼠源抗原结构）。

---

## 🧰 扩展模块 1：人-鼠抗原表位一致性验证模块

**目的**：确保人源抗原和鼠源抗原的结构配对不仅是同源蛋白，而是在抗体结合表位区域具有结构或序列一致性。

**输入**：

* 阶段四中对鼠合的结构配对（每个 pair 有 human_ag_chain.cif 和 mouse_ag_chain_aligned.cif）
* 阶段一元数据（含抗体链信息）

**处理逻辑**：

1. **提取人源抗原的表位残域**：

   * 加载处理合结构（人源抗原 + 抗体），识别所有抗原链上与抗体链距离小于 5Å 的残域；
   * 使用 Bio.PDB 或 MDAnalysis 脚本计算原子-原子间距离。

2. **映射这些残域到鼠源抗原**：

   * 先进行序列比对（pairwise2 或 Clustal）建立残域索引映射；
   * 在鼠源结构中提取对应残域坐标。

3. **评估结构一致性**：

   * 计算对应残域间 RMSD；
   * 记录匹配残域个数、序列 identity 百分比、结构 RMSD（Cα）值。

4. **结果写入 CSV**：为每个人-鼠结构对添加字段：

   * 表位_identity（%）
   * 表位_RMSD（Å）
   * 是否为高一致性（自定义判定标准，如 RMSD < 1.5Å 且 identity > 80%）

**使用库**：Bio.PDB、Bio.Align、numpy、pandas

**输出**：增强版 CSV 文件（如 human_mouse_pairs.csv）中添加上述字段。

---

## 🤐 扩展模块 2：SAbDab 鼠源抗原复合物识别模块

**目的**：在尝试 AlphaFold DB 前，优先查询 SAbDab 是否已经存在该人源抗原的鼠源同源结构。

**输入**：

* SAbDab 元数据表（含抗原名称、PDB ID、抗原链、抗原物种）
* 当前目标人源抗原（PDB ID + chain）对应的抗原名称或 UniProt ID

**处理逻辑**：

1. **查找同名或同源抗原条目**：

   * 在 SAbDab 元数据中过滤所有 `抗原名称相同` 或 `UniProt ID 同源` 且 `抗原物种为 mouse` 的结构条目；

2. **检查这些鼠源条目是否为抗体-抗原复合物**：

   * 若是，则提取其鼠源抗原链并作为替代模型使用。

3. **作为鼠源结构输入**：

   * 若找到该结构，直接用于结构对齐；
   * 如果多个结构候选，则优先选分辨率高、结构完整、与人源抗原链类型一致者。

**优点**：

* 使用真实复合物结构更可靠；
* 若该复合物中已包含抗体，则可直接从该结构中识别鼠源抗原表位，用于跨种对比。

**输出**：替代 AlphaFold 模型的 `mouse_ag_chain.cif`，并记录来源标注为 `SAbDab`（而非 `AlphaFold`）

**使用库**：pandas、Bio.PDB、requests（如需下载结构）

**增量机制**：记录已查找过的抗原名称组合，避免重复搜索。将结构来源写入配对元数据表中。

---

其余五个阶段维持原结构，在结构获取和验证阶段中引入上述两个模块作为补充判断逻辑。这样整个流程既遵守“不可AI建模”的限制，又保证人-鼠结构对在抗体表位层面具有真正的一致性。
