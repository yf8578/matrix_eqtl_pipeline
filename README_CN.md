# MatrixEQTL 分析流程

> [English Version (英文版)](README.md)

这一流程整合了 PLINK (`process_vcf_plink.py`) 和 MatrixEQTL (`run_matrix_eqtl.py`)，提供了一套从 VCF 到最终 eQTL 结果的完整、稳健的分析方案。

## 🌟 主要特性

* **极速处理**：使用 PLINK 进行 VCF 过滤和格式转换，比纯 Python 快数十倍。
* **内存优化**：支持大文件分块读取 (`chunk-size`)。
* **完整协变量**：自动计算基因型 PCA (人群分层) 和 PEER 因子 (隐藏混杂因素)。
* **双模式**：支持基于位置的 Cis/Trans 分析 (Mode A) 和 无位置信息的全基因组扫描 (Mode B)。
* **自动质控**：内置正态性检查 (`check_normality.py`) 和协变量相关性检查 (`check_covariates.py`)。

---

## 🛠️ 环境配置 (Environment Setup)

为确保流程顺利运行，建议使用 Conda 构建独立环境：

### 1. 创建基础环境

```bash
# 创建包含 Python 3.9 和 R 4.x 的环境
conda create -n eqtl_env python=3.9 r-base=4.3
conda activate eqtl_env

# 安装 PLINK (需添加 bioconda channel)
conda install -c bioconda plink
```

### 2. 安装 Python 依赖

```bash
# 安装分析与绘图所需的包
pip install pandas scikit-learn rpy2 matplotlib seaborn
```

### 3. 安装 R 依赖

```bash
# 安装 MatrixEQTL
R -e 'install.packages("MatrixEQTL", repos="http://cran.us.r-project.org")'

# 安装 PEER (推荐使用仓库内自带的源码包)
# 方法 A: 源码安装 (Linux/Mac)
R CMD INSTALL R_peer_source_1.3.tgz

# 方法 B: Conda 方式 (作为备选)
conda install -c bioconda r-peer
```

---

## 🚀 快速开始 (Quick Start)

### 1. 设置环境变量

将以下代码保存为脚本或直接在终端运行，请**务必修改**为你自己的数据路径：

```bash
# === 1. 设置路径 (请修改这里) ===
export WORK_DIR=$(pwd)
export CODE_DIR=${WORK_DIR}/matrix_eqtl_pipeline/code

# 输入文件路径
export VCF_FILE="/path/to/your/genotypes.vcf.gz"
export EXP_FILE="/path/to/your/expression_matrix.tsv"
# 基因位置文件 (Mode A 必需)
export GENE_LOC_FILE="/path/to/your/gene_locations.tsv"

# 创建输出目录
mkdir -p ${WORK_DIR}/output
```

### 2. 第一步：VCF 处理 (使用 PLINK)

这一步会清洗 VCF，保留与表达量矩阵重叠的样本，并生成用于 PCA 的二进制文件和用于 MatrixEQTL 的矩阵。

```bash
echo ">>> Step 1: Processing VCF..."
python3 ${CODE_DIR}/process_vcf_plink.py \
    --vcf ${VCF_FILE} \
    --expression ${EXP_FILE} \
    --out-dir ${WORK_DIR}/output/step1_plink \
    --maf 0.05 \
    --plink-cmd plink \
    --snps-only \
    --threads 4
```

### 3. 第二步：生成协变量 (PCA + PEER)

```bash
echo ">>> Step 2: Generating Covariates..."

# 2.1 表达量正态化 (IQN)
python3 ${CODE_DIR}/run_iqn.py \
    -i ${EXP_FILE} \
    -o ${WORK_DIR}/output/step2_covariates/expression.qnorm

# 2.2 基因型 PCA (PLINK 计算前 3 个主成分)
plink \
    --bfile ${WORK_DIR}/output/step1_plink/plink_temp \
    --pca 3 \
    --threads 4 \
    --out ${WORK_DIR}/output/step2_covariates/plink_pca

# 格式化 PCA 结果
python3 ${CODE_DIR}/format_plink_pca.py \
    -i ${WORK_DIR}/output/step2_covariates/plink_pca.eigenvec \
    -o ${WORK_DIR}/output/step2_covariates/genotype.pcs \
    -n 3

# 2.3 PEER 隐因子挖掘 (比如计算 15 个因子)
# 如果有已知协变量(如性别)，可加上 -c known_covariates.txt
python3 ${CODE_DIR}/run_peer.py \
    -i ${WORK_DIR}/output/step2_covariates/expression.qnorm \
    -n 15 \
    -o ${WORK_DIR}/output/step2_covariates/peer_factors.tsv

# 2.4 合并所有协变量
python3 ${CODE_DIR}/combine_covariates.py \
    -p ${WORK_DIR}/output/step2_covariates/genotype.pcs \
    -f ${WORK_DIR}/output/step2_covariates/peer_factors.tsv \
    -o ${WORK_DIR}/output/step2_covariates/final_covariates.txt
```

### 4. 第三步：运行 MatrixEQTL 分析

#### 选项 A: 标准 eQTL (区分 Cis/Trans)

适用于你有基因位置信息，且关注 Cis 调控的情况。需要 `gene_positions.txt` (表头: geneid, chr, left, right)。

```bash
echo ">>> Step 3: Running Mode A (Cis/Trans)..."

# (可选) 修正基因位置文件表头
# sed '1s/features/geneid/' ${GENE_LOC_FILE} > ${WORK_DIR}/output/step1_plink/gene_positions.txt

python3 ${CODE_DIR}/run_matrix_eqtl.py \
    --genotype-matrix ${WORK_DIR}/output/step1_plink/genotype.matrix \
    --genotype-positions ${WORK_DIR}/output/step1_plink/snp_positions.txt \
    --gene-expression-matrix ${WORK_DIR}/output/step2_covariates/expression.qnorm \
    --gene-positions ${WORK_DIR}/output/step1_plink/gene_positions.txt \
    --covariates ${WORK_DIR}/output/step2_covariates/final_covariates.txt \
    --output-file ${WORK_DIR}/output/final_cis/cis_results.txt \
    --trans-output-file ${WORK_DIR}/output/final_cis/trans_results.txt \
    --p-value 0.05 \
    --trans-p-value 1e-5 \
    --cis-distance 1000000 \
    --chunk-size 2000
```

#### 选项 B: 全局分析 (无位置信息 / mQTL)

适用于没有具体染色体位置，或单纯做全基因组关联扫描。

```bash
echo ">>> Step 4: Running Mode B (Location-Agnostic)..."
python3 ${CODE_DIR}/run_matrix_eqtl_noloc.py \
    --genotype-matrix ${WORK_DIR}/output/step1_plink/genotype.matrix \
    --expression-matrix ${WORK_DIR}/output/step2_covariates/expression.qnorm \
    --covariates ${WORK_DIR}/output/step2_covariates/final_covariates.txt \
    --output-file ${WORK_DIR}/output/final_noloc/global_results.txt \
    --p-value 1e-5 \
    --no-fdr \
    --chunk-size 5000
```

### 5. 第四步 (可选): 使用 QTLtools 进行分析

如果你追求更严谨的 FDR 校正（Permutation）或需要进行条件分析，可以使用本流程集成的 QTLtools 模块。

> **官方网站**: [https://qtltools.github.io/qtltools/](https://qtltools.github.io/qtltools/)

我们提供了简化的 Python 封装脚本，帮你自动处理格式转换和命令构建。

#### 5.1 准备数据 (生成 BED 格式)

**场景 A: 有基因位置 (标准 eQTL)**

```bash
# 自动合并表达量和位置文件，并进行 bgzip/tabix 压缩索引
python3 ${CODE_DIR}/qtltools_prep.py \
    --expression ${WORK_DIR}/output/step2_covariates/expression.qnorm \
    --positions ${GENE_LOC_FILE} \
    --out ${WORK_DIR}/output/qtltools_input/expression.bed
```

**场景 B: 无基因位置 (脂质/代谢物)**

```bash
# 生成伪位置 (默认 chr1:1000 起)
python3 ${CODE_DIR}/qtltools_prep.py \
    --expression ${EXP_FILE} \
    --dummy-pos \
    --out ${WORK_DIR}/output/qtltools_input/lipid_dummy.bed
```

*注意：还需要手动准备协变量文件（如果你有的话），可以直接复用 Step 2 的结果：*

```bash
cp ${WORK_DIR}/output/step2_covariates/final_covariates.txt ${WORK_DIR}/output/qtltools_input/covariates.txt
```

#### 5.2 运行分析 (Permutation / Nominal / Trans)

使用统一的 `run_qtltools.py` 脚本：

```bash
# 例子：运行 Cis Permutation (推荐用于 eQTL)
python3 ${CODE_DIR}/run_qtltools.py \
    --vcf ${VCF_FILE} \
    --bed ${WORK_DIR}/output/qtltools_input/expression.bed.gz \
    --cov ${WORK_DIR}/output/qtltools_input/covariates.txt \
    --out ${WORK_DIR}/output/qtltools_results/permutations.txt \
    --mode permute \
    --permutations 1000

# 例子：运行 Trans 分析 (推荐用于代谢物/无位置性状)
python3 ${CODE_DIR}/run_qtltools.py \
    --vcf ${VCF_FILE} \
    --bed ${WORK_DIR}/output/qtltools_input/lipid_dummy.bed.gz \
    --cov ${WORK_DIR}/output/qtltools_input/covariates.txt \
    --out ${WORK_DIR}/output/qtltools_results/trans_results.txt \
    --mode trans
```

---
---

## 📜 脚本功能详解

### 核心脚本

* `process_vcf_plink.py` (Step 1)
  * **功能**: 一站式 VCF 处理。
  * **特点**: 自动处理样本重叠，保留 SNP，多线程加速。
  * **参数**: `--vcf`, `--expression`, `--maf`(默认0.05), `--threads`.

* `run_matrix_eqtl.py` (Step 3 - Mode A)
  * **功能**: 运行带位置信息的标准 Cis/Trans eQTL 分析。
  * **关键参数**: `--cis-distance` (定义Cis范围), `--p-value` (Cis阈值), `--trans-p-value` (Trans阈值).

* `run_matrix_eqtl_noloc.py` (Step 3 - Mode B)
  * **功能**: 运行无位置信息的 MatrixEQTL (ANOVA 模型或全局扫描)。
  * **特点**: 适合 mQTL 或未知位置的转录本。

### 辅助脚本

* `run_iqn.py`
  * **功能**: 表达量矩阵的正态化 (Inverse Quantile Normalization)。
  * **更新**: 运行后会自动生成一张 `.check.png` 图片，用于检查数据是否符合正态分布。

* `run_peer.py`
  * **功能**: 计算 PEER 隐因子。
  * **更新**: 支持 `-c/--covariates` 参数，可以把已知协变量（如性别、年龄）先传入模型进行校正。

* `format_plink_pca.py`
  * **功能**: 将 PLINK 的 `.eigenvec` 文件转置为 MatrixEQTL 需要的格式。

* `check_covariates.py` (新增质控)
  * **功能**: 计算你的 **已知协变量** 与 **PEER因子** / **表达量PCs** 之间的相关性。
  * **用途**: 帮助你判断哪些协变量是冗余的，或者验证批次效应是否被捕捉到。

---

## 📂 附录：文件格式示例

### 1. MatrixEQTL

**输入：基因型/表达量矩阵 (`.txt` / `.tsv`)**

```tsv
id      sample1 sample2 sample3 ...
SNP_01  0       1       2       ...
SNP_02  1       1       0       ...
```

**输入：位置文件 (`.txt`)**

```tsv
geneid  chr     left    right
Gene_A  1       100     200
Gene_B  X       500     600
```

*(注意：MatrixEQTL 不关心表头名称，但顺序必须是 ID, Chr, Pos1, Pos2)*

**输出：结果文件 (`results.txt`)**

```tsv
SNP     gene    beta    t-stat  p-value FDR
SNP_01  Gene_A  0.5     2.3     1e-5    0.01
```

### 2. QTLtools

**输入：表型文件 (`.bed`)**
*(由 `qtltools_prep.py` 自动生成)*

```tsv
#Chr    start   end     pid     gid     strand  sample1 sample2 ...
1       999     1000    Gene_A  Gene_A  .       1.2     -0.5    ...
1       1999    2000    Gene_B  Gene_B  .       0.8     2.1     ...
```

**输出：Permutation (`.txt`)**

```tsv
phe_id  n_var   mle_shape1  mle_shape2  dummy   rel_beta    adj_emp_pval ...
Gene_A  500     ...         ...         ...     ...         1.2e-4
```

*(重点关注最后一列 `adj_emp_pval`，即校正后的 P 值)*

**输出：Nominal (`.txt`)**

```tsv
phe_id  var_id  dist    pval    slope   is_best
Gene_A  SNP_01  123     1e-5    0.5     1
```
