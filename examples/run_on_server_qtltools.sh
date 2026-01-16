#!/bin/bash
# =============================================================================
# QTLtools Analysis Master Template (run_on_server_qtltools.sh)
# =============================================================================
# This script runs the full pipeline adapted for QTLtools.
# Steps:
#   Step 1: Sample Intersection
#   Step 2: Expression Processing (Norm -> BED.gz)
#   Step 3: Genotype Processing (VCF.gz)
#   Step 4: Covariates (PCA + PEER -> TXT)
#   Step 5: QTLtools Analysis (Cis Permutation / Nominal)
# =============================================================================

# --- Configuration ---
export PROJ_DIR='/data/work/QTLpipeline'
export CODE_DIR="/data/work/QTLpipeline/matrix_eqtl_pipeline"
export OUTPUT_DIR="$PROJ_DIR/results_qtltools"

# Input Data
VCF_FILE="/data/work/SNP_data/control/only_control_20240613_v2.vcf.gz"
EXPRESSION_FILE="/data/work/xQTL/pQTL/01_data/control/49_control_protein_expression_matrix.tsv"
GTF_FILE="/data/work/Ref/gencode.v46.annotation.gtf" # Optional if using fetch-pos
KNOWN_COV_FILE="/data/work/new_QTL/output/step2_covariates/20260107_sample_info_for_zyf.csv"
QTLTOOLS_BIN="QTLtools"  # Path to QTLtools binary

# Parameters
THREADS=4
PLINK_BIN="plink2"   # Or "plink" for v1.9
PLINK_VERSION=2      # 1 or 2

echo ">>> Project Dir: $PROJ_DIR"
echo ">>> Output Dir:  $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# Step 1: Sample Intersection
# =============================================================================
echo "[Step 1] Finding Common Samples..."
python3 "$CODE_DIR/src/matrix_eqtl/utils/io.py" \
    --vcf "$VCF_FILE" \
    --expression "$EXPRESSION_FILE" \
    --out "$OUTPUT_DIR/samples.txt"

# =============================================================================
# Step 2: Expression Processing (Norm -> BED.gz)
# =============================================================================
echo "[Step 2] Processing Expression (BED Format)..."

# 2a: Filter & Normalize (Output: .qnorm)
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/expression.py" \
    --expression "$EXPRESSION_FILE" \
    --samples "$OUTPUT_DIR/samples.txt" \
    --out-dir "$OUTPUT_DIR/preprocessing" \
    --no-tpm-filter \
    --max-missing 0.1 \
    # --no-norm (Optional: remove comment to skip INT)

# 2b: Convert to QTLtools BED (Merge with Positions + Compress)
# Uses the .qnorm file from 2a
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/locations.py" \
    --expression "$OUTPUT_DIR/preprocessing/expression.qnorm" \
    --out "$OUTPUT_DIR/preprocessing/expression.bed" \
    --format qtltools \
    --fetch-pos --id-type uniprot_gn_id,entrezgene_id,external_gene_name,ensembl_gene_id
    # Note: This produces expression.bed.gz and expression.bed.gz.tbi

# =============================================================================
# Step 3: Genotype Processing (VCF.gz)
# =============================================================================
echo "[Step 3] Processing Genotypes (Filtered VCF)..."
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/genotype.py" \
    --vcf "$VCF_FILE" \
    --samples "$OUTPUT_DIR/samples.txt" \
    --out-dir "$OUTPUT_DIR/preprocessing" \
    --tool qtltools \
    --maf 0.05 \
    --geno 0.05 \
    --hwe 1e-6 \
    --threads "$THREADS" \
    --plink-bin "$PLINK_BIN" \
    --plink-version "$PLINK_VERSION"
    # Output: genotypes.filtered.vcf.gz

# =============================================================================
# Step 4: Covariates (PCA + PEER + Combine)
# =============================================================================
echo "[Step 4] Generating Covariates..."

# 4a: Genotype PCA (Using MatrixEQTL-style preprocessing for speed/simplicity)
# Note: Genotype PCA is same regardless of downstream tool
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/genotype_pca.py" \
    --genotype "$OUTPUT_DIR/preprocessing/genotypes" \
    --out-dir "$OUTPUT_DIR/covariates" \
    --pca-n 3 \
    --threads "$THREADS" \
    --plink-bin "$PLINK_BIN" \
    --plink-version "$PLINK_VERSION"

# 4b: PEER Factors (Same logic)
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/expression_peer.py" \
    --expression "$OUTPUT_DIR/preprocessing/expression.qnorm" \
    --out-dir "$OUTPUT_DIR/covariates" \
    --peer-n 5

# 4c: Combine
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/combine_covariates.py" \
    --pca "$OUTPUT_DIR/covariates/genotype_pcs.txt" \
    --peer "$OUTPUT_DIR/covariates/peer_factors.tsv" \
    --known "$KNOWN_COV_FILE" \
    --out-dir "$OUTPUT_DIR/covariates"

# Output is final_covariates.txt (Tab separated). QTLtools accepts this.

# =============================================================================
# Step 5: QTLtools Analysis (Official Conditional Analysis Workflow)
# =============================================================================
echo "[Step 5] Running QTLtools Analysis..."
echo "  Workflow: Permutation -> FDR Thresholds -> Conditional Analysis"

BED_FILE="$OUTPUT_DIR/preprocessing/expression.bed.gz"
VCF_FILE_QC="$OUTPUT_DIR/preprocessing/genotypes.filtered.vcf.gz"
COV_FILE="$OUTPUT_DIR/covariates/final_covariates.txt"

# -----------------------------------------------------------------------------
# 5a: Permutation Pass
# Calculate adjusted P-values for every phenotype
# -----------------------------------------------------------------------------
echo "  [5a] Running Cis Permutations (1000 iter)..."
"$QTLTOOLS_BIN" cis \
    --vcf "$VCF_FILE_QC" \
    --bed "$BED_FILE" \
    --cov "$COV_FILE" \
    --permute 10000 \
    --normal \
    --std-err \
    --out "$OUTPUT_DIR/permutations_all.txt"

# -----------------------------------------------------------------------------
# 5b: Calculate FDR Thresholds
# Use the included R script to determine P-value threshold for FDR < 0.05
# -----------------------------------------------------------------------------
echo "  [5b] Calculating FDR Thresholds (FDR=0.05)..."
# Compress permutations file for R script
gzip -f "$OUTPUT_DIR/permutations_all.txt"

# Locate runFDR_cis.R (Assuming it's in the same dir as QTLtools or provided path)
# If not found, you need to download it from QTLtools website/repo
# Here we assume it is in $CODE_DIR/r_src or similar. 
# Providing a placeholder path - ADJUST THIS PATH to where runFDR_cis.R actually is
FDR_SCRIPT="$CODE_DIR/r_src/qtltools_runFDR_cis.R" 

# if [ ! -f "$FDR_SCRIPT" ]; then
#     echo "Warning: $FDR_SCRIPT not found. Cannot run FDR step."
#     echo "Please download qtltools_runFDR_cis.R from https://github.com/qtltools/qtltools/tree/master/scripts"
#     exit 1
# fi

Rscript "$FDR_SCRIPT" \
    "$OUTPUT_DIR/permutations_all.txt.gz" \
    0.05 \
    "$OUTPUT_DIR/permutations_all"

# This produces: $OUTPUT_DIR/permutations_all.thresholds.txt

# -----------------------------------------------------------------------------
# 5c: Conditional Analysis
# Use the thresholds to find all independent signals
# -----------------------------------------------------------------------------
echo "  [5c] Running Conditional Analysis..."
"$QTLTOOLS_BIN" cis \
    --vcf "$VCF_FILE_QC" \
    --bed "$BED_FILE" \
    --cov "$COV_FILE" \
    --mapping "$OUTPUT_DIR/permutations_all.thresholds.txt" \
    --normal \
    --std-err \
    --out "$OUTPUT_DIR/conditional_all.txt"

# -----------------------------------------------------------------------------
# 5d: Extract Top Variants
# Filter for rank 0 (first independent signal) or significant conditional hits
# -----------------------------------------------------------------------------
echo "  [5d] Extracting Significant Results..."
# Extracting primary signals (rank 0) or all independent signals?
# The doc says: awk '{ if ($21 == 1) print $0 }' for backward p-value significance
awk '{ if ($21 == 1) print $0 }' "$OUTPUT_DIR/conditional_all.txt" > "$OUTPUT_DIR/sig_results_conditional.txt"

echo ">>> QTLtools Pipeline Finished."
echo "  Results: $OUTPUT_DIR/sig_results_conditional.txt"
