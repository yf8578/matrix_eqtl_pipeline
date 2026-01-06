# MatrixEQTL & QTLtools 分析流程

> [English Version (英文版)](README.md)

这是一个整合了 **MatrixEQTL** (快速矩阵运算) 和 **QTLtools** (Permutation/FDR 校正) 的综合 eQTL/pQTL 分析流程。它旨在提供一站式的解决方案，从 VCF 基因型和表达量矩阵开始，自动处理数据质控、格式转换、协变量计算，并生成最终的关联分析结果。

---

## 🌟 主要特性

* **双引擎支持**：
  * **MatrixEQTL**：基于 R 的矩阵运算，速度极快，支持 Cis/Trans (Mode A) 和无位置扫描 (Mode B)。
  * **QTLtools**：基于 C++ 的标准工具，支持 Permutation 分析 (FDR 校正) 和条件分析。
* **批量自动化**：提供 `template.sh` 和 CSV 配置系统，可一次性生成大量独立任务脚本。
* **流程解耦**：MatrixEQTL 和 QTLtools 流程完全独立，互不依赖，输出目录隔离。
* **自动协变量**：内置 PLINK PCA (人群分层) 和 PEER (隐因子) 自动计算流程。
* **极速预处理**：使用 PLINK 进行 VCF 过滤，比纯 Python 处理快数十倍。

---

## 🛠️ 环境配置

建议使用 Conda 构建独立环境：

```bash
# 1. 创建环境
conda create -n eqtl_env python=3.9 r-base=4.3
conda activate eqtl_env

# 2. 安装核心工具
conda install -c bioconda plink qtltools tabix

# 3. 安装 Python 依赖
pip install pandas scikit-learn rpy2 matplotlib seaborn

# 4. 安装 R 包
R -e 'install.packages("MatrixEQTL", repos="http://cran.us.r-project.org")'
# 安装 PEER (需自行编译或通过 conda r-peer)
```

---

## 🚀 批量分析系统 (推荐)

我们提供了一套批量脚本生成系统，适合同时分析多个组织、多个性状或不同参数组合。

### 1. 准备 CSV 配置文件

创建一个 CSV 文件 (例如 `batch_input.csv`)，每一行代表一个独立的分析任务：

| id | work_dir | vcf | expression | location | peer | cis_dist |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| blood_cis | /data/work/blood | /data/vcf/chr1.vcf.gz | /data/exp/blood.tsv | /data/genes.loc | 15 | 1000000 |
| brain_cis | /data/work/brain | /data/vcf/chr1.vcf.gz | /data/exp/brain.tsv | /data/genes.loc | 10 | 1000000 |

* `id`: 任务ID，用于命名脚本。
* `location`: 基因位置文件。如果为空，则只会生成 MatrixEQTL 无位置模式脚本，不会生成 QTLtools 脚本。
* `peer`: PEER 因子数量 (默认 15)。

### 2. 生成运行脚本

运行 `template.sh`：

```bash
bash template.sh batch_input.csv
```

程序会在 `generated_scripts/` 目录下为每个任务生成**两套独立脚本**：

1. **`run_job_{id}_matrix.sh`**: MatrixEQTL 分析流程
    * 输出目录: `{work_dir}/output_matrix/`
2. **`run_job_{id}_qtltools.sh`**: QTLtools 分析流程
    * 输出目录: `{work_dir}/output_qtltools/`
    * *仅当提供了 `location` 时生成*

### 3. 提交任务

可以直接运行生成的脚本，或者使用 `run_all_jobs.sh` 串行运行所有任务：

```bash
# 运行单个任务
bash generated_scripts/run_job_blood_cis_qtltools.sh

# 或 提交到集群 (视集群环境而定)
# qsub generated_scripts/run_job_blood_cis_qtltools.sh
```

---

## 📖 手动分步指南 (Manual Guide)

如果您想手动运行或了解内部流程，可以参考以下步骤。我们提供了独立的脚本来演示完整流程。

### 方式 A: MatrixEQTL 流程

主要用于快速扫描或 Cis/Trans 共分析。

1. **VCF 预处理**: `process_vcf_plink.py` 使用 PLINK 提取重叠样本。
2. **生成协变量**: `run_iqn.py` (正态化) -> PLINK PCA -> `run_peer.py` -> 合并。
3. **运行分析**: `run_matrix_eqtl.py` (含位置) 或 `run_matrix_eqtl_noloc.py` (无位置)。

### 方式 B: QTLtools 流程 (Cis Permutation)

主要用于需要严格 FDR 校正的 Cis-eQTL 研究。

**核心步骤命令示例：**

```bash
# 1. 准备 BED 文件 (使用 qtltools_prep.py 辅助)
# 该脚本自动处理：
#   - 合并表达量和位置文件
#   - 处理表头 (保证 #chr start end pid gid strand 格式)
#   - 排序基因组坐标
#   - bgzip 压缩和 tabix 索引
python3 code/qtltools_prep.py \
    --expression expression.qnorm \
    --positions gene_locations.tsv \
    --out expression.bed

# 2. 准备协变量 (可以直接复用之前步骤生成的，或自行准备)
cp final_covariates.txt covariates.txt

# 3. 运行 QTLtools (直接调用命令)
QTLtools cis \
    --vcf genotypes.vcf.gz \
    --bed expression.bed.gz \
    --cov covariates.txt \
    --out permutations.txt \
    --permute 1000 \
    --window 1000000 \
    --normal
```

*注意：我们也提供了一个独立的 `run_qtltools_standalone.sh` 脚本，包含了从 VCF 处理到最终 QTLtools 分析的全过程，您可以直接复制修改使用。*

---

## 📂 脚本说明

* `code/create_batch_scripts.py`: 批量脚本生成器核心逻辑。
* `code/qtltools_prep.py`: 专门用于生成 QTLtools 标准 BED 格式的工具。
* `code/run_qtltools.py`: (旧版) Python 封装的 QTLtools 运行器，现在推荐在脚本中直接写 `QTLtools` 命令以便于定制。
* `code/mxeqtl.R`: MatrixEQTL 的 R 核心代码 (已修复 helper function 缺失问题)。
