# 目标
根据项目的要求产生对应的 yaml 文件
在每一个簇里面，选择 reference 的 antibody 作为参考；设计 CDR 区域让其对所有 antigen 都具有不错的免疫原性。
这里是要生成符合要求的 yaml 文件以便使用。

## 实现细节 (Implementation Details)

### 1. 数据来源 (Data Sources)
*   **簇目录 (Cluster Directories)**: `data/cleaned_split/{cluster_id}/`。每个 `{cluster_id}` 文件夹代表一个簇，其中包含多个 `{epitope_id}/` 子文件夹。
*   **簇成员文件 (Cluster Member Files)**: `data/cleaned_split/{cluster_id}/{epitope_id}/` 目录下包含 `antibody.cif` 和 `antigen.cif` 文件。
*   **簇元数据 (Cluster Metadata)**: `data/cluster/epitope_cluster_summary.csv` 包含簇中每个 `epitope_id` 的详细信息。
*   **SAbDab 总结文件 (SAbDab Summary)**: `data/raw_data/meta/sabdab_summary_all.tsv` 用于查找 PDB 条目的物种信息，以便选择参考抗体。
*   **表位信息 (Epitope Information)**: `data/extract_epitope/pockets/{epitope_id}.npz` 包含每个表位的残基级别信息。

### 2. 参考抗体选择逻辑 (Reference Antibody Selection Logic)

对于每个簇：
1.  收集该簇中所有 `epitope_id` 对应的 PDB ID。
2.  从 `sabdab_summary_all.tsv` 中查找这些 PDB 的物种 (`organism` 或 `heavy_species`/`light_species` 字段)。
3.  **优先级**:
    *   首先选择**人源 (Homo sapiens)** 抗体。
    *   如果没有人源抗体，则选择**鼠源 (Mus musculus)** 抗体。
    *   如果两者都没有，则选择簇中列出的**第一个** PDB 作为参考。
4.  确定参考抗体的 `epitope_id` (例如 `6nf2_Q`)。

### 3. YAML 文件生成逻辑 (YAML File Generation Logic)

每个簇生成一个 YAML 文件，保存到 `data/yamls/boltzgen_like/multi_antigen/cluster_{cluster_id}.yaml`。

**YAML 结构详情 (Detailed YAML Structure)**:

```yaml
multi_antigen:
  antibody:
    file:
      path: "data/cleaned_split/{cluster_id}/{reference_epitope_id}/antibody.cif"
      include:
        # 从 reference_epitope_id/antibody.cif 文件头部的 _entity_poly.pdbx_strand_id 中解析重链和轻链的 chain ID
        - chain: {id: "U"} # 例如：重链的 Chain ID
        - chain: {id: "V"} # 例如：轻链的 Chain ID
    chains: ["U", "V"] # 与 include 中的 Chain ID 保持一致
  antigens:
    # 遍历该簇目录 (data/cleaned_split/{cluster_id}/) 下的所有 {member_epitope_id} 文件夹
    - name: "{member_epitope_id}" # 例如："6nf2_Q"
      file:
        path: "data/cleaned_split/{cluster_id}/{member_epitope_id}/antigen.cif"
        include:
          # 从 member_epitope_id/antigen.cif 文件头部的 _entity_poly.pdbx_strand_id 中解析抗原链的 chain ID
          - chain: {id: "Q"} # 例如：抗原链的 Chain ID
      epitope:
        # 从 data/extract_epitope/pockets/{member_epitope_id}.npz 中解析表位残基信息
        - chain: "Q" # 抗原链的 Chain ID
          residue_ids: [12, 13, 14, 15, 16] # 表位残基的序列号列表 (residue sequence numbers)
  parallel_config:
    aggregation_method: "kabsch_mean"
    antigen_weights: null # 保持为 null
  # constraints 字段暂时不添加，根据用户指示，目前抗体设计用不到此字段
```

### 4. 关键步骤和考量 (Key Steps and Considerations)

*   **路径处理**: YAML 文件中的所有文件路径都将使用相对于项目根目录的路径 (`data/...`)。
*   **链 ID 解析**: 需要编写代码来解析 `antibody.cif` 和 `antigen.cif` 文件，以提取其 `_entity_poly.pdbx_strand_id` 中的实际链 ID。这是因为链 ID 可能不是标准的 'H', 'L', 'A' 等。
*   **表位残基提取**: 将加载 `.npz` 文件 (例如 `data/extract_epitope/pockets/{epitope_id}.npz`)，从中提取 `chain_ids` 和 `res_seq_nums` 数组，并整理成 `{chain: [residue_ids]}` 的格式。

---
```python
    def valid_yaml_content(self):
        """Create valid multi_antigen YAML content."""
        return {
            "multi_antigen": {
                "antibody": {
                    "file": {
                        "path": "antibody.cif",
                        "include": [
                            {"chain": {"id": "H"}},
                            {"chain": {"id": "L"}},
                        ],
                    },
                    "chains": ["H", "L"],
                },
                "antigens": [
                    {
                        "name": "Antigen1",
                        "file": {
                            "path": "antigen1.cif",
                            "include": [{"chain": {"id": "A"}}],
                        },
                    },
                    {
                        "name": "Antigen2",
                        "file": {
                            "path": "antigen2.cif",
                            "include": [{"chain": {"id": "B"}}],
                        },
                    },
                ],
                "parallel_config": {
                    "aggregation_method": "kabsch_mean",
                    "antigen_weights": None,
                },
            },
            "constraints": [
                {"bond": {"atom1": ["H", 22, "SG"], "atom2": ["H", 92, "SG"]}}
            ],
        }
```