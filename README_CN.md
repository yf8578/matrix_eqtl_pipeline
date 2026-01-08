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

## 📖 手动分步指南 (Step-by-Step Guide)

如果您想手动运行每一步，或者想了解脚本背后的具体操作，请仔细阅读以下流程。

### 第一阶段：通用数据准备 (Step 1-2)
所有软件都需要完成以下两步：VCF处理和协变量计算。

#### Step 1: VCF 预处理
**目的**：从庞大的 VCF 中提取出与表达量表重叠的样本，并转换为二进制格式以加速后续计算。
```bash
python3 code/process_vcf_plink.py \
    --vcf "/path/to/genotypes.vcf.gz" \
    --expression "/path/to/expression.tsv" \
    --out-dir "output/step1_plink" \
    --maf 0.05 \
    --threads 4 \
    --snps-only
# 输出: output/step1_plink/plink_temp.{bed,bim,fam} (提取后的SNP数据)
```

#### Step 2: 计算协变量 (Covariates)
**目的**：校正混杂因素（如批次效应、人群结构），提高分析的统计效力。无需手动计算，直接运行以下脚本：
```bash
# 2.1 表达量正态化 (IQN)
python3 code/run_iqn.py -i "/path/to/expression.tsv" -o "output/step2_covariates/expression.qnorm"

# 2.2 基因型 PCA (人群结构)
plink --bfile "output/step1_plink/plink_temp" --pca 3 --out "output/step2_covariates/plink_pca"
python3 code/format_plink_pca.py -i "output/step2_covariates/plink_pca.eigenvec" -o "output/step2_covariates/genotype.pcs" -n 3

# 2.3 PEER 隐因子 (批次效应)
python3 code/run_peer.py -i "output/step2_covariates/expression.qnorm" -n 15 -o "output/step2_covariates/peer_factors.tsv"

# 2.4 合并最终协变量文件
python3 code/combine_covariates.py \
    -p "output/step2_covariates/genotype.pcs" \
    -f "output/step2_covariates/peer_factors.tsv" \
    -o "output/step2_covariates/final_covariates.txt"
```

---

### 第二阶段：关联分析 (Step 3) - 三选一

根据您的需求，选择以下三种引擎之一进行后续分析。

#### 🅰️ 选项 A: 使用 MatrixEQTL (CPU, 极快)
适用于：初看结果、Trans-eQTL 筛选。

1. **准备位置文件** (如有需要)：
```bash
# 自动抓取位置
python3 code/prep_genes.py --expression "output/step2_covariates/expression.qnorm" --fetch-pos --id-type ensembl_gene_id --out "output/gene_locations.tsv" --format matrix
```

2. **运行 MatrixEQTL**：
```bash
python3 code/run_matrix_eqtl.py \
    --genotype-matrix "output/step1_plink/genotype.matrix" \
    --gene-expression-matrix "output/step2_covariates/expression.qnorm" \
    --covariates "output/step2_covariates/final_covariates.txt" \
    --gene-positions "output/gene_locations.tsv" \
    --output-file "output/matrix_results.txt" \
    --cis-distance 1000000
```
*结果文件*：`output/matrix_results.txt` (包含 p-value, FDR, beta 等)。

---

#### 🅱️ 选项 B: 使用 QTLtools (CPU, 严谨)
适用于：正式发表文章 (Cis-eQTL Permutation Analysis)。

1. **准备 BED 输入文件**：
**必须**使用 `prep_genes.py` 生成包含位置信息的压缩 BED 文件：
```bash
python3 code/prep_genes.py \
    --expression "output/step2_covariates/expression.qnorm" \
    --fetch-pos --id-type ensembl_gene_id \
    --vcf "path/to/vcf" \
    --out "output/qtltools_input/expression.bed"
# 输出: output/qtltools_input/expression.bed.gz (和 .tbi 索引)
```

2. **运行 QTLtools Permutation**：
```bash
QTLtools cis \
    --vcf "path/to/vcf" \
    --bed "output/qtltools_input/expression.bed.gz" \
    --cov "output/step2_covariates/final_covariates.txt" \
    --out "output/qtltools_results/permutations.txt" \
    --window 1000000 --permute 1000 --normal
```
*结果文件*：`output/qtltools_results/permutations.txt`。主要关注 `adj_empirical_pval` 列。

---

#### 🚀 选项 C: 使用 TensorQTL (GPU, 极速且严谨)
适用于：有 GPU 的服务器，追求极致速度与精确性。

1. **准备输入文件**：
*   同 QTLtools，需要 `expression.bed.gz` (由 prep_genes.py 生成)。
*   需要 PLINK 文件 `output/step1_plink/plink_temp` (由 Step 1 生成)。

2. **运行 TensorQTL**：
使用我们提供的 Python 包装器：
```bash
python3 code/run_tensorqtl.py \
    --plink "output/step1_plink/plink_temp" \
    --phenotypes "output/qtltools_input/expression.bed.gz" \
    --covariates "output/step2_covariates/final_covariates.txt" \
    --mode cis \
    --out "output/tensorqtl_results/cis"
```
*结果文件*：`output/tensorqtl_results/cis.cis_perm.txt.gz`。格式与 FastQTL 兼容。

---

### 3. `create_batch_scripts.py` (批量生成)

**功能**：读取 CSV 配置文件，批量生成 MatrixEQTL / QTLtools / TensorQTL 的运行脚本。

```bash
# 1. 准备 BED 文件 (使用 prep_genes.py 辅助)
# 该脚本自动处理：
#   - 合并表达量和位置文件
#   - 处理表头 (保证 #chr start end pid gid strand 格式)
#   - 排序基因组坐标
#   - bgzip 压缩和 tabix 索引
python3 code/prep_genes.py \
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
* `code/prep_genes.py`: 专门用于生成 QTLtools 标准 BED 格式或 MatrixEQTL 位置文件的全能工具。
* `code/run_qtltools.py`: (旧版) Python 封装的 QTLtools 运行器，现在推荐在脚本中直接写 `QTLtools` 命令以便于定制。
* `code/mxeqtl.R`: MatrixEQTL 的 R 核心代码 (已修复 helper function 缺失问题)。

---

### 4. 软件对比与选择指南 (Which one to choose?)

| 特性 | MatrixEQTL | QTLtools | TensorQTL |
| :--- | :--- | :--- | :--- |
| **语言/平台** | R (CPU) | C++ (CPU) | Python + PyTorch (GPU/CPU) |
| **主要优势** | **极其简单好用**，适合初学者。Trans-eQTL 极快。 | **学术界金标准**。提供 Permutation 和 Beta 近似，结果最严谨。 | **速度之王 (GPU)**。算法与 QTLtools 一致，但快 100-300 倍。 |
| **核心算法** | 线性回归 (Linear Regression), ANOVA | 线性回归 + 置换检验 (Permutation) | 线性回归 + 置换检验 (Permutation) |
| **适用场景** | 小样本，或只做 Trans-eQTL 筛选。 | 中型样本，或需要投稿高水平期刊 (Reviewer 认可度高)。 | **大样本 (>500)**，或有 GPU 服务器。Biobank 首选。 |
| **缺点** | 默认校正 (FDR) 对 LD 结构处理不如 Permutation 严谨。 | 无 GPU 加速，大样本跑 Cis-Permutation 极慢。 | 需要配置 CUDA 环境，依赖 PyTorch。 |

#### 📊 哪个更可靠？

* **可靠性 (Reliability)**: **QTLtools 和 TensorQTL 并列第一**。
  * QTLtools 是目前的行业标准 (Gold Standard)。
  * TensorQTL 是 QTLtools 的 GPU 复刻版，数学模型完全一致 (GTEx V8 官方推荐)。
* **MatrixEQTL** 也很可靠，但它的 P 值校正逻辑相对简单，在 Cis-eQTL 分析中不如前两者严谨。

#### ✅ 推荐方案

1. **首选 TensorQTL**: 只要服务器有显卡，它就是最佳选择（又快又准）。
2. **次选 QTLtools**: 如果没显卡，但要发文章，选它最稳妥。
3. **备选 MatrixEQTL**: 如果只想快速看个大概，或者做 Trans 分析。
