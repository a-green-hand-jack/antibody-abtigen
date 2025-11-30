# 03_generate.py 脚本说明

该脚本的主要功能是根据去重后的 CSV 数据文件，批量生成用于后续处理（如 Boltz 推理或特征提取）的 YAML 配置文件。

## 核心功能

1.  **读取数据**：
    *   从指定的 CSV 文件（通常是经过 `02_decuplication.py` 处理后的数据）读取抗体-抗原复合物信息。

2.  **构建 YAML 结构**：
    *   遍历每一行数据，提取重链（H_chain）、轻链（L_chain，如果存在）和抗原链（antigen_chain）的信息。
    *   构建符合特定格式（类似 Boltz 输入格式）的字典结构：
        ```yaml
        version: 1
        sequences:
          - protein:
              id: <Chain ID>
              sequence: <Sequence>
          ...
        ```

3.  **批量保存**：
    *   将每个条目保存为单独的 `.yaml` 文件。
    *   文件名以 CSV 中的 `file_name` 字段命名。

## 主要函数

*   `convert_single_entry_to_yaml(entry_value, version=1)`: 将单行 DataFrame 数据转换为 YAML 对应的字典结构。
*   `convert_to_yaml(entry_key, entry_yaml_content, save_folder)`: 使用 `ruamel.yaml` 库将字典写入文件，并保持特定的缩进格式。
*   `main(args)`: 主流程，负责目录创建、数据遍历和文件生成。

## 使用参数

*   `--processed_data_dir_after_decuplication`: 去重后的 CSV 文件路径（默认：`./data/summary-decuplication-distance_threshold_9.csv`）。
*   `--output_dir`: YAML 文件输出目录（默认：`./data/yaml_for_data_after_decuplication-distance_threshold_9`）。

## 注意事项
*   脚本会检查输出目录是否存在，如果存在会**先删除再重新创建**，请注意备份数据。
