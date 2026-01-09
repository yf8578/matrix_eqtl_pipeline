# MatrixEQTL 统一分析流程 (Unified MatrixEQTL Pipeline)

> [English Version (英文版)](README.md)

这是一个整合了 **MatrixEQTL**、**QTLtools** 和 **TensorQTL** 的统一 eQTL 分析流程。本文档重点介绍 **MatrixEQTL** 的详细手动分析步骤。

---

## 🔍 手动分步执行指南 (Manual Step-by-Step Guide - MatrixEQTL)

根据最新的代码更新，以下是针对 MatrixEQTL 分析的详细每一步操作指令。这些指令包含了对大文件处理、自定义过滤和多线程的支持。

### Step 0: 准备位置文件 (Prepare Location Files)

MatrixEQTL 需要两个位置文件：基因位置 (`gene_locations.txt`) 和 SNP 位置 (`snp_positions.txt`)。
SNP 位置会在 **Step 3** 自动生成，因此这一步只需生成基因位置。

**命令**:

```bash
python3 src/matrix_eqtl/preprocessing/locations.py \
    --expression data/expression.tsv \
    --out data/gene_locations.txt \
    --format matrixeqtl \
    --fetch-pos \
    --id-type ensembl_gene_id
```

* `--fetch-pos`: 自动联网从 BioMart 下载基因位置 (推荐)。
* `--gtf your.gtf`: 或者指定本地 GTF 文件解析。
* `--id-type`: 基因 ID 类型 (默认 `ensembl_gene_id`, 可选 `external_gene_name`, `uniprot_gn_id` 等)。

### Step 1: 样本对齐 (Sample Intersection)

找出 VCF 和表达量的共有样本。如果无共有样本，脚本会报错并打印各自的前5个ID供检查。

**命令**:

```bash
python3 src/matrix_eqtl/utils/io.py \
    --vcf data/genotypes.vcf.gz \
    --expression data/expression.tsv \
    --out samples.txt
```

### Step 2: 表达量预处理 (Expression Processing)

对表达量矩阵进行过滤和标准化。

**命令**:

```bash
python3 src/matrix_eqtl/preprocessing/expression.py \
    --expression data/expression.tsv \
    --samples samples.txt \
    --out-dir results/preprocessing \
    --no-tpm-filter \
    --max-missing 0.1 \
    --no-norm
```

**参数详解**:

* `--no-tpm-filter`: **(新增)** 跳过基于数值的过滤（如 `TPM > 0.1`）。适用于已经处理过的数值（如蛋白数据）或只需要基于缺失值过滤的情况。
* `--max-missing 0.1`: **(新增)** 剔除缺失值 (NaN) 比例超过 10% 的基因。
* `--min-tpm 0.1`: (如果不使用 `--no-tpm-filter`) 定义最小表达阈值。
* `--min-samples-rel 0.2`: (如果不使用 `--no-tpm-filter`) 定义需满足阈值的最小样本比例。
* `--no-norm`: (可选) 跳过 Inverse Normal Transformation (INT) 标准化。

### Step 3: 基因型预处理 (Genotype Processing)

使用 PLINK 进行质控 (QC) 并转换为 MatrixEQTL 格式。此步骤已针对千万级变异位点优化，支持流式处理以防内存溢出。

**命令**:

```bash
python3 src/matrix_eqtl/preprocessing/genotype.py \
    --vcf data/genotypes.vcf.gz \
    --samples samples.txt \
    --out-dir results/preprocessing \
    --tool matrixeqtl \
    --maf 0.01 \
    --geno 0.05 \
    --hwe 1e-6 \
    --threads 4
```

**参数详解**:

* `--maf 0.01`: 最小等位基因频率阈值。
* `--geno 0.05`: SNP 缺失率阈值。
* `--hwe 1e-6`: 哈代-温伯格平衡 p 值阈值。
* `--threads 4`: **(新增)** PLINK 使用的线程数，加速 QC 和转码步骤。
* 此步骤会自动生成 `genotype.matrix` (转置后的剂量矩阵) 和 `snp_positions.txt`。

### Step 4: 协变量生成 (Covariate Generation)

运行 PCA (Genotype) 和 PEER (Expression) 生成协变量。

**命令**:

```bash
python3 src/matrix_eqtl/preprocessing/covariates.py \
    --genotype results/preprocessing/genotypes \
    --expression results/preprocessing/expression.qnorm \
    --known data/known_covariates.txt \
    --out-dir results/covariates \
    --pca-n 3 \
    --peer-n 5 \
    --threads 4
```

**参数详解**:

* `--known`: (可选) 已知协变量文件。
* `--pca-n 3`: 基因型主成分数量 (由 PLINK 计算)。
* `--peer-n 5`: 表达量 PEER 因子数量 (由 R `peer` 包计算)。
* `--threads 4`: **(新增)** PLINK PCA 使用的线程数。
* **注意**: 如果 `--genotype` 文件不存在，脚本会报错并列出目录内容以供调试。

### Step 5: 执行 MatrixEQTL 分析 (Run Analysis)

最后运行 MatrixEQTL 核心分析。

**命令**:

```bash
python3 src/matrix_eqtl/analysis/run_matrix_eqtl.py \
    --genotype-matrix results/preprocessing/genotype.matrix \
    --gene-expression-matrix results/preprocessing/expression.qnorm \
    --covariates results/covariates/final_covariates.txt \
    --gene-positions data/gene_locations.txt \
    --genotype-positions results/preprocessing/snp_positions.txt \
    --output-file results/cis_results.txt \
    --cis-distance 1000000 \
    --p-value 1e-5
```

**参数详解**:

* `--cis-distance`: Cis-eQTL 的窗口大小 (bp)，通常为 1000000 (1MB)。
* `--p-value`: 输出结果的 P 值阈值。

---

## 🛠️ 其他工具用法 (QTLtools & TensorQTL)

## 🔍 进阶分析与参数说明 (Advanced Analysis)

### 1. Cis 与 Trans 单独阈值控制

MatrixEQTL 支持同时分析 Cis 和 Trans eQTL，并设置不同的 P 值阈值。

* `--cis-distance`: Cis 窗口大小 (默认为 1e6 即 1MB)。
* `--p-value`: **CIS** 分析的 P 值阈值 (如果不提供 `--trans-p-value`，则默认也用于 Trans)。
* `--trans-p-value`: **TRANS** 分析的 P 值阈值 (默认为 0，即不保留 Trans 结果)。
* `--trans-output-file`: Trans 结果的输出文件路径。

**示例：同时保留 Cis (1e-5) 和 Trans (1e-8) 结果**

```bash
python3 src/matrix_eqtl/analysis/run_matrix_eqtl.py \
    ... \
    --p-value 1e-5 \
    --output-file results/cis_results.txt \
    --trans-p-value 1e-8 \
    --trans-output-file results/trans_results.txt
```

### 2. FDR 校正 (False Discovery Rate)

MatrixEQTL 默认会计算 FDR (Benjamini-Hochberg) 并输出在结果的最后一列。

* 默认行为：输出结果包含 FDR 列。
* `--no-fdr`: **(新增参数)** 如果数据量极大导致内存不足，或不需要 FDR，使用此参数关闭 FDR 计算。

    ```bash
    python3 src/matrix_eqtl/analysis/run_matrix_eqtl.py ... --no-fdr
    ```

### 3. 无位置文件分析 (Location-Free / ANOVA Mode)

如果您的数据无法区分 Cis/Trans (例如：没有位置信息，或想进行全基因组两两关联)，可以不提供 `--gene-positions` 和 `--genotype-positions` 参数。

**行为**:

* 所有的 Gene-SNP 对都会被测试 (Cis-distance 参数被忽略)。
* 分析模式等同于全 Trans 分析。

**命令示例**:

```bash
python3 src/matrix_eqtl/analysis/run_matrix_eqtl.py \
    --genotype-matrix results/preprocessing/genotype.matrix \
    --gene-expression-matrix results/preprocessing/expression.qnorm \
    --covariates results/covariates/final_covariates.txt \
    --output-file results/all_pairs_results.txt \
    --p-value 1e-10  # 全基因组分析建议使用更严格的阈值
```

---

## 🛠️ 其他工具用法 (QTLtools & TensorQTL)

(待补充 - 流程类似，主要区别在于 Step 0 的位置文件格式和 Step 3 的输出格式)
