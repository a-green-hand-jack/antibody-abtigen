将 PDB/CIF 文件转化为 `boltzgen` (基于 Boltz-1 的生成/预测工具) 的 YAML 配置文件，核心在于**数据清洗**和**意图定义**。

CIF 文件是\*\*“结果”**（包含坐标、实验元数据），而 YAML 文件是**“指令”\*\*（告诉模型生成什么、限制什么）。

要把 CIF 转换为 YAML，你需要完成以下三个核心步骤的映射。下面我将针对**抗原-抗体复合物**（特别是 scFv/Fab）提供具体的转换逻辑和 Python 代码示例。

### 1\. 核心转换逻辑

#### A. 序列 (Sequences) - **最重要**

  * **CIF 源**: `_entity_poly.pdbx_seq_one_letter_code` (必须读取这个，**不能**只读取有坐标的原子，因为 Loop 和 Linker 往往在坐标里是缺失的)。
  * **YAML 目标**: `sequences` 列表。
  * **注意**: 对于 scFv，你要确保提取的是**包含 Linker 的完整长序列**。

#### B. 二硫键 (Disulfide Bonds)

  * **CIF 源**: `_struct_conn` 类别，其中 `conn_type_id` 为 `disulf`。
  * **YAML 目标**: `bonds` 列表。
  * **难点**: CIF 中的二硫键通常使用的是 PDB 编号（如 Cys H22），你需要将其转换为序列的**索引位置 (0-indexed)**。

#### C. Loop / 缺失区域 / 设计区域

  * **策略**: 如果你是做**预测 (Prediction)**，只需提供完整序列，模型会自动补全 Loop。如果你是做**设计 (Design)**（比如想重新生成 scFv 的 CDR 区），你需要定义 `design_id` 或掩码。

-----

### 2\. Python 自动化转换脚本 (使用 `gemmi`)

我为你写了一个 Python 脚本，使用 `gemmi` 库（处理 CIF 最专业的库）来提取信息并生成 YAML 格式。

你需要先安装库：`pip install gemmi pyyaml`

```python
import gemmi
import yaml

def cif_to_boltzgen_yaml(cif_path, output_yaml_path):
    # 读取 CIF 文件
    doc = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()

    # 1. 提取序列信息 (Entities)
    # -------------------------------------------------
    sequences = []
    entity_id_to_seq = {}
    
    # 遍历 _entity_poly 类别获取完整序列
    for row in block.find('_entity_poly.', ['entity_id', 'pdbx_seq_one_letter_code']):
        entity_id = row['entity_id']
        # 移除换行符，获取纯序列
        raw_seq = row['pdbx_seq_one_letter_code'].replace('\n', '')
        # 将3字母代码转换为单字母 (CIF有时存的是带括号的非标准格式，这里简化处理，假设是标准单字母)
        # 如果序列中有非标准氨基酸，可能需要额外处理
        seq = raw_seq 
        entity_id_to_seq[entity_id] = seq
        
        sequences.append({
            "protein": {
                "id": f"chain_{entity_id}",  # 这里的 ID 可以自定义
                "sequence": seq
            }
        })

    # 2. 提取二硫键 (Disulfide Bonds)
    # -------------------------------------------------
    bonds = []
    
    # 获取结构连接信息
    struct_conn = block.find('_struct_conn.', 
                             ['conn_type_id', 
                              'ptnr1_auth_seq_id', 'ptnr1_auth_comp_id', # 原子1
                              'ptnr2_auth_seq_id', 'ptnr2_auth_comp_id'  # 原子2
                             ])
    
    # 注意：这里有一个复杂点。CIF记录的是 "Residue Number" (如 22, 104)，
    # 但 YAML 通常需要序列中的 "Index" (第几个氨基酸，从0或1开始)。
    # 严格的转换需要将 PDB ID 映射回 Sequence Index。
    # 下面的代码是一个简化的逻辑，演示如何提取二硫键对。
    
    for row in struct_conn:
        if row['conn_type_id'] == 'disulf':
            # 这里我们只能提取到类似 "23" 和 "88" 这样的残基编号
            # 在实际工程中，你需要建立一个 PDB_ID -> Seq_Index 的映射表
            bond_info = {
                "atom1_idx": row['ptnr1_auth_seq_id'], # 需转换为 index
                "atom2_idx": row['ptnr2_auth_seq_id'], # 需转换为 index
                "type": "disulfide"
            }
            bonds.append(bond_info)

    # 3. 构建最终字典
    # -------------------------------------------------
    boltz_data = {
        "sequences": sequences,
        # 如果有二硫键，取消注释下面这行
        # "bonds": bonds 
    }
    
    # 4. 写入 YAML
    # -------------------------------------------------
    with open(output_yaml_path, 'w') as f:
        yaml.dump(boltz_data, f, sort_keys=False, default_flow_style=False)
    
    print(f"转换完成！已保存至 {output_yaml_path}")
    print(f"提取了 {len(sequences)} 条链。")

# 使用示例 (替换为你的文件名)
# cif_to_boltzgen_yaml("9BSB.cif", "9BSB_design.yaml")
```

### 3\. 手动调整指南 (针对抗原-抗体复合物)

生成的 YAML 只是一个基础骨架，针对你的 **scFv 抗体设计** 需求，你需要手动检查和修改以下几个关键部分：

#### A. 序列部分 (`sequences`)

如果是 scFv（如 9BSB）：

```yaml
sequences:
  - protein:
      id: A
      # 确保这里包含：VH序列 + Linker(GGGGS...) + VL序列
      # 即使 CIF 结构里 Linker 缺失，这里必须写上，否则模型会把它们预测成断开的两个蛋白。
      sequence: "QVQLQ...[VH]...GGGGSGGGGS...[Linker]...DIQMT...[VL]"
```

#### B. 键的定义 (`bonds`)

抗体内部的二硫键对稳定性至关重要。虽然 AlphaFold/Boltz 能够隐式学习到，但显式指定通常更好。
你需要查找抗体的保守位点（通常是 Cys23 和 Cys104，卡巴编号）。

```yaml
bonds:
  - [0, 22, 96] # 示例：第一条链(0)，第22号残基 和 第96号残基 形成二硫键
```

*注意：Boltz 的索引方式可能不同（全局索引 vs 链内索引），请参照具体文档。*

#### C. 设计 vs 固定 (`constraints` / `generation`)

这是 `boltzgen` 的核心。如果你想保留抗原，设计抗体，你需要设置约束。

假设 Chain A 是抗原 (Antigen)，Chain B 是 scFv (Antibody)。

**场景 1：全自动预测 (Structure Prediction)**
不需要额外设置，直接跑。它会自动预测 Linker 的形态。

**场景 2：固定抗原，优化抗体 CDR (Binder Design)**
你需要告诉工具哪些部分是“死的”，哪些是“活的”。

```yaml
# 这是一个概念性的示例，具体语法参考 boltzgen 文档
constraints:
  - type: "fixed_structure"
    chain: "A"  # 锁死抗原结构
    path: "antigen.pdb" # 提供抗原的参考结构

design:
  - target: "B" # 针对 scFv
    # 如果你想让 AI 重新设计 Linker，你可以把 Linker 序列变成 Mask
    # 或者指定需要采样的区域
```

### 总结

1.  **Linker 处理**：在 YAML `sequence` 里必须**手动补全** Linker 的氨基酸序列。不要让它空着，否则 scFv 会断成两截。
2.  **Loop 处理**：对于 CDR Loop 或其他柔性 Loop，**不要**在序列中删除它们。只需提供完整序列，AI 模型（Boltz-1/AlphaFold3）的一个主要功能就是预测这些看不见的 Loop 的构象。
3.  **二硫键**：如果你的 scFv 容易解折叠，显式地在 YAML 中定义二硫键会有帮助。
4.  **从小分子入手**：如果你的复合物里有小分子（如 LSD），别忘了在 `sequences` 列表里加上 `ligand` 或 `smiles` 定义，因为抗体结合构象往往依赖于配体的存在。
