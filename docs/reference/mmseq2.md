# MMseqs2 安装与配置指南 (GPU 加速版)

本文档详细说明如何使用 Mamba 创建一个独立的环境来管理 MMseqs2 (GPU 版本)，并将预编译的二进制文件部署到指定路径，以便跨项目共享。

> **注意**：官方的 Bioconda/Conda-forge 源中的 `mmseqs2` 通常**只包含 CPU 版本**。我们将手动下载 GPU 版本的二进制文件并放置在 Mamba 环境的 `bin` 目录下。

## 1. 环境准备 (使用 Mamba)

我们将创建一个名为 `mmseqs` 的环境，并将其放置在 `/ibex/user/wuj0c/envs/mmseqs` 目录下。

```bash
# 1. 创建一个新的环境（指定路径）
mamba create -p /ibex/user/wuj0c/envs/mmseqs python=3.10 -y

# 2. 激活环境
conda activate /ibex/user/wuj0c/envs/mmseqs
# 或者
source activate /ibex/user/wuj0c/envs/mmseqs

# 3. 安装 CUDA Toolkit
# 这里的 11.8 是一个兼容性较好的版本
mamba install -c conda-forge cudatoolkit=11.8 -y
```

## 2. 安装 MMseqs2 (GPU 版本)

我们将直接下载官方预编译的 GPU 二进制文件，并将其放置在环境的 `bin` 目录下，这样在激活环境后可以直接调用 `mmseqs`，或者通过绝对路径调用。

1.  **进入环境的 bin 目录**：
    ```bash
    cd /ibex/user/wuj0c/envs/mmseqs/bin
    ```

2.  **下载并解压预编译包**：
    访问 [MMseqs2 GitHub Releases](https://github.com/soedinglab/MMseqs2/releases) 确认下载链接。以下以 `linux-gpu` 版本为例：

    ```bash
    # 下载最新的 GPU 版本 (Linux)
    wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz
    
    # 解压
    tar xvfz mmseqs-linux-gpu.tar.gz
    
    # 将解压出的 mmseqs 二进制文件移动到当前 bin 目录
    # 注意：解压后通常会在当前目录下生成一个 mmseqs 文件夹
    mv mmseqs/bin/mmseqs .
    
    # 赋予执行权限
    chmod +x mmseqs
    
    # 清理
    rm mmseqs-linux-gpu.tar.gz
    rm -rf mmseqs # 删除解压出来的其余文件夹
    ```

3.  **验证安装**：
    ```bash
    ./mmseqs --help
    ```

## 3. 验证 GPU 支持

安装完成后，可以通过以下命令检查 MMseqs2 是否正确识别了 GPU。

```bash
# 运行 help 命令
/ibex/user/wuj0c/envs/mmseqs/bin/mmseqs --help

# 或者尝试运行一个简单的 search 命令并观察显存使用
# /ibex/user/wuj0c/envs/mmseqs/bin/mmseqs search ... --gpu 1
```

## 4. 在脚本中使用

在运行 `05_cluster.py` 时，请通过 `--mmseqs_path` 参数指定到刚刚安装的绝对路径。

```bash
python scripts/reference/05_cluster.py \
    --processed_csv_path data/summary-decuplication-distance_threshold_9.csv \
    --out_dir data \
    --mmseqs_path /ibex/user/wuj0c/envs/mmseqs/bin/mmseqs \
    --seq_identity_threshold 0.5
```

### 脚本参数说明
脚本 `05_cluster.py` 默认并未开启 `--gpu 1`。如果需要强制开启 GPU 加速，请确保在调用 `subprocess.run` 时添加了 `--gpu 1` 参数，或者确认 `easy-cluster` 在检测到 GPU 时是否自动启用。建议显式修改脚本以确保使用 GPU。

## 5. 实际部署记录 (2025-11-30)

以下是在 `/ibex/user/wuj0c/envs/mmseqs` 实际部署时的操作记录，特别注意解决了解压后文件夹名称冲突的问题。

```bash
# 1. 创建环境
mamba create -p /ibex/user/wuj0c/envs/mmseqs python=3.10 -y

# 2. 安装依赖
mamba install -p /ibex/user/wuj0c/envs/mmseqs -c conda-forge cudatoolkit=11.8 -y

# 3. 下载与安装 MMseqs2
cd /ibex/user/wuj0c/envs/mmseqs/bin
wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz
tar xvfz mmseqs-linux-gpu.tar.gz

# 4. 处理文件冲突并移动二进制文件
# 解压出的文件夹也叫 mmseqs，不能直接 mv mmseqs/bin/mmseqs . (因为当前目录下已有同名目录)
cp mmseqs/bin/mmseqs ./mmseqs_bin   # 先复制并重命名
chmod +x mmseqs_bin                 # 赋予执行权限
rm -rf mmseqs                       # 删除解压出的文件夹
mv mmseqs_bin mmseqs                # 改回原名

# 5. 清理
rm mmseqs-linux-gpu.tar.gz
```

