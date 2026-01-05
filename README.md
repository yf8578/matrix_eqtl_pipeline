# MatrixEQTL Analysis Pipeline

> [中文文档 (Chinese Version)](README_CN.md)

This pipeline integrates **PLINK** (`process_vcf_plink.py`) and **MatrixEQTL** (`run_matrix_eqtl.py`) to provide a robust, high-performance solution for eQTL analysis, scaling from VCF processing to final association testing.

## 🌟 Key Features

* **High Performance**: Uses PLINK for VCF filtering and conversion, significantly faster than Python-based parsing.
* **Memory Optimized**: Supports chunked processing for large datasets.
* **Comprehensive Covariates**: Automatically generates Genotype PCA (for population structure) and PEER factors (for hidden confounders).
* **Dual Modes**: Supports **Mode A** (Cis/Trans analysis with location data) and **Mode B** (Location-agnostic / Global scan).
* **Automated QC**: Includes tools for normality checks (`check_normality.py`) and covariate correlation checks (`check_covariates.py`).

---

## 🚀 Quick Start

### 0. Prerequisites

Ensure the following are installed and in your PATH:

* **PLINK** (v1.9 or v2.0)
* **Python 3** (with `pandas`, `sklearn`, `rpy2`)
* **R** (with `MatrixEQTL`, `peer` packages)

### 1. Set Up Variables

Copy and paste this block into your terminal (modify paths for your data):

```bash
# === 1. Set Paths (MODIFY THESE) ===
export WORK_DIR=$(pwd)
export CODE_DIR=${WORK_DIR}/matrix_eqtl_pipeline/code

# Input Files
export VCF_FILE="/path/to/your/genotypes.vcf.gz"
export EXP_FILE="/path/to/your/expression_matrix.tsv"
# Gene Locations (Required for Mode A)
export GENE_LOC_FILE="/path/to/your/gene_locations.tsv"

# Create Output Directory
mkdir -p ${WORK_DIR}/output
```

### 2. Step 1: VCF Processing (via PLINK)

Cleans VCF, filters for sample overlap, retains only SNPs, and generates formats for PCA and MatrixEQTL.

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

### 3. Step 2: Generate Covariates (PCA + PEER)

```bash
echo ">>> Step 2: Generating Covariates..."

# 2.1 Expression Normalization (IQN)
python3 ${CODE_DIR}/run_iqn.py \
    -i ${EXP_FILE} \
    -o ${WORK_DIR}/output/step2_covariates/expression.qnorm

# 2.2 Genotype PCA (Calculate top 3 PCs using PLINK)
plink \
    --bfile ${WORK_DIR}/output/step1_plink/plink_temp \
    --pca 3 \
    --threads 4 \
    --out ${WORK_DIR}/output/step2_covariates/plink_pca

# Format PCA output for MatrixEQTL
python3 ${CODE_DIR}/format_plink_pca.py \
    -i ${WORK_DIR}/output/step2_covariates/plink_pca.eigenvec \
    -o ${WORK_DIR}/output/step2_covariates/genotype.pcs \
    -n 3

# 2.3 PEER Factor Estimation (e.g., 15 factors)
# Optional: Add -c known_covariates.txt to include known factors (Sex, Age)
python3 ${CODE_DIR}/run_peer.py \
    -i ${WORK_DIR}/output/step2_covariates/expression.qnorm \
    -n 15 \
    -o ${WORK_DIR}/output/step2_covariates/peer_factors.tsv

# 2.4 Combine All Covariates
python3 ${CODE_DIR}/combine_covariates.py \
    -p ${WORK_DIR}/output/step2_covariates/genotype.pcs \
    -f ${WORK_DIR}/output/step2_covariates/peer_factors.tsv \
    -o ${WORK_DIR}/output/step2_covariates/final_covariates.txt
```

### 4. Step 3: Run Analysis

#### Option A: Standard eQTL (Cis/Trans)

Use this if you have gene/SNP positions and want to distinguish **Cis** (local) vs **Trans** (distant) effects.

```bash
echo ">>> Step 3: Running Mode A (Cis/Trans)..."

# (Optional) Fix gene header: MatrixEQTL requires 'geneid'
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

#### Option B: Location-Agnostic / Global Scan

Use this for mQTL or if physical distance is not relevant.

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

---

## 📜 Detailed Script Reference

### Core Scripts

* `process_vcf_plink.py` (Step 1)
  * **Purpose**: All-in-one VCF processing using PLINK.
  * **Features**: Handles overlapping, SNP filtering, multithreading (`--threads`).

* `run_matrix_eqtl.py` (Step 3 - Mode A)
  * **Purpose**: Wrapper for MatrixEQTL with location information (Cis/Trans).
  * **Key Args**: `--cis-distance` (default 1Mb), `--p-value` (cis threshold), `--trans-p-value`.

* `run_matrix_eqtl_noloc.py` (Step 3 - Mode B)
  * **Purpose**: Location-agnostic MatrixEQTL (ANOVA model or global scan).
  * **Features**: Optimized for speed and large datasets.

### Helper Utilities

* `run_iqn.py`
  * **Purpose**: Inverse Quantile Normalization of expression data.
  * **Update**: Automatically plots a normality check (`.check.png`).

* `run_peer.py`
  * **Purpose**: Calculates PEER hidden factors.
  * **Update**: Supports `-c/--covariates` to include known covariates (e.g. Sex) in the model.

* `format_plink_pca.py`
  * **Purpose**: Converts PLINK `.eigenvec` to MatrixEQTL format.

* `check_covariates.py` (New QA Tool)
  * **Purpose**: Checks correlations between **Known Covariates** and **PEER Factors** / **Expression PCs**.
  * **Usage**: Run this to validate if your covariates are relevant or redundant.

* `combine_covariates.py`
  * Combines genotype PCs, PEER factors, and additional covariates into one file.

*(Legacy scripts like `pc_covariates.py` are preserved but PLINK PCA is recommended)*
