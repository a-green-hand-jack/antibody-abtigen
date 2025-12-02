我们现在要调整聚类的思路，之前的操作可以参考[对齐结构](./11_align_structures.md)和[使用表位聚类](./10_cluster_epitopes.md)

我们现在要修改聚类的思路，核心思想是“从数据里，找出抗体相同的记录，然后比对抗原的结构和序列相似性，找出哪些你需要的抗原”

具体的讲我们使用 [summary文件](../../data/raw_data/meta/sabdab_summary_all.tsv)确定哪些抗原是对应一个抗体，把哪些具有一个抗体的抗原-抗体复合物作为一个簇。

然后继续使用 “对齐结构”根据之前已经处理好的[表位数据](data/extract_epitope)得到新的对齐的 cif 文件。

然后再根据之前已经有的清理 cif 文件的[脚本](scripts/13_clean_and_split.py)去得到干净的 cif。

我们这个脚本本质上和 [根据表位聚类](docs/scripts/10_cluster_by_same_antibody.md)是同级别的，目的是一样的，但是操作不一。
