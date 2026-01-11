#!/bin/bash
# =============================================================================
# MatrixEQTL Analysis Master Template (run_on_server.sh)
# =============================================================================
# This script runs the full pipeline from raw data to final eQTL results.
# It includes:
#   Step 0: Location Files Preparation
#   Step 1: Sample Intersection (VCF vs Expression)
#   Step 2: Expression Processing (Filtering + Normalization)
#   Step 3: Genotype Processing (PLINK QC + Formatting)
#   Step 4: Covariate Generation (PCA + PEER + Combine)
#   Step 5: eQTL Analysis (Standard Cis/Trans & Location-Free)
# =============================================================================

# --- Configuration (Edit these paths) ---
export PROJ_DIR=$(pwd)
export CODE_DIR="$PROJ_DIR"  # Path to the source code (src/)
export OUTPUT_DIR="$PROJ_DIR/results_control"

# Input Data
VCF_FILE="/data/work/SNP_data/control/only_control_20240613_v2.vcf.gz"
EXPRESSION_FILE="/data/work/xQTL/pQTL/01_data/control/49_control_protein_expression_matrix.tsv"
GTF_FILE="/data/work/Ref/gencode.v46.annotation.gtf" # Optional if using fetch-pos
KNOWN_COV_FILE="/data/work/new_QTL/output/step2_covariates/20260107_sample_info_for_zyf.csv"

# Parameters
THREADS=4
CIS_DIST=1000000     # 1MB
CIS_PVAL=1e-5
TRANS_PVAL=1e-5      # Set >0 to output trans results. If 0, no trans output.
ALL_PAIRS_PVAL=1e-10 # For location-free mode (usually stricter)

echo ">>> Project Dir: $PROJ_DIR"
echo ">>> Output Dir:  $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# Step 0: Prepare Location Files
# =============================================================================
echo "-----------------------------------------------------------------------------"
echo "[Step 0] Checking Location Files..."
if [ ! -f "$OUTPUT_DIR/gene_locations.txt" ]; then
    echo "Generating gene locations..."
    python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/locations.py" \
        --expression "$EXPRESSION_FILE" \
        --out "$OUTPUT_DIR/gene_locations.txt" \
        --format matrixeqtl \
        --fetch-pos --id-type ensembl_gene_id
else
    echo "Gene locations already exist at $OUTPUT_DIR/gene_locations.txt"
fi

# =============================================================================
# Step 1: Sample Intersection
# =============================================================================
echo "-----------------------------------------------------------------------------"
echo "[Step 1] Finding Common Samples..."
python3 "$CODE_DIR/src/matrix_eqtl/utils/io.py" \
    --vcf "$VCF_FILE" \
    --expression "$EXPRESSION_FILE" \
    --out "$OUTPUT_DIR/samples.txt"
    
# =============================================================================
# Step 2: Expression Processing
# =============================================================================
echo "-----------------------------------------------------------------------------"
echo "[Step 2] Processing Expression Data..."
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/expression.py" \
    --expression "$EXPRESSION_FILE" \
    --samples "$OUTPUT_DIR/samples.txt" \
    --out-dir "$OUTPUT_DIR/preprocessing" \
    --no-tpm-filter \
    --max-missing 0.1 \
    # --no-norm
    # Remove --no-norm if you want Inverse Normal Transformation (INT)

# =============================================================================
# Step 3: Genotype Processing
# =============================================================================
echo "-----------------------------------------------------------------------------"
echo "[Step 3] Processing Genotypes (PLINK QC)..."
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/genotype.py" \
    --vcf "$VCF_FILE" \
    --samples "$OUTPUT_DIR/samples.txt" \
    --out-dir "$OUTPUT_DIR/preprocessing" \
    --tool matrixeqtl \
    --maf 0.05 \
    --geno 0.05 \
    --hwe 1e-6 \
    --threads "$THREADS"

# =============================================================================
# Step 4: Covariate Generation (Split Steps)
# =============================================================================
echo "-----------------------------------------------------------------------------"
echo "[Step 4] Generating Covariates..."

# 4a: Genotype PCA
echo "  [4a] Running Genotype PCA..."
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/genotype_pca.py" \
    --genotype "$OUTPUT_DIR/preprocessing/genotypes" \
    --out-dir "$OUTPUT_DIR/covariates" \
    --pca-n 3 \
    --threads "$THREADS"

# 4b: Expression PEER
echo "  [4b] Calculating PEER Factors..."
# Check if R peer is available, otherwise skip or handle error
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/expression_peer.py" \
    --expression "$OUTPUT_DIR/preprocessing/expression.qnorm" \
    --out-dir "$OUTPUT_DIR/covariates" \
    --peer-n 5

# 4c: Combine
echo "  [4c] Combining Covariates..."
python3 "$CODE_DIR/src/matrix_eqtl/preprocessing/combine_covariates.py" \
    --pca "$OUTPUT_DIR/covariates/genotype_pcs.txt" \
    --peer "$OUTPUT_DIR/covariates/peer_factors.tsv" \
    --known "$KNOWN_COV_FILE" \
    --out-dir "$OUTPUT_DIR/covariates"

# Check correlation plot
if [ -f "$OUTPUT_DIR/covariates/covariates_correlation.pdf" ]; then
    echo "  >> Correlation plot generated: $OUTPUT_DIR/covariates/covariates_correlation.pdf"
fi

# =============================================================================
# Step 5: Run Analysis (Two Modes)
# =============================================================================

# --- Mode A: Standard Cis + Trans Analysis (With Location Files) ---
echo "-----------------------------------------------------------------------------"
echo "[Step 5a] Running Standard Cis + Trans Analysis..."

# NOTE: We output both Cis and Trans results here.
# Cis distance is set by --cis-distance.
# Trans results are enabled by setting --trans-p-value > 0.

python3 "$CODE_DIR/src/matrix_eqtl/analysis/run_matrix_eqtl.py" \
    --genotype-matrix "$OUTPUT_DIR/preprocessing/genotype.matrix" \
    --gene-expression-matrix "$OUTPUT_DIR/preprocessing/expression.qnorm" \
    --covariates "$OUTPUT_DIR/covariates/final_covariates.txt" \
    --gene-positions "$OUTPUT_DIR/gene_locations.txt" \
    --genotype-positions "$OUTPUT_DIR/preprocessing/snp_positions.txt" \
    --output-file "$OUTPUT_DIR/cis_results.txt" \
    --cis-distance "$CIS_DIST" \
    --p-value "$CIS_PVAL" \
    --trans-output-file "$OUTPUT_DIR/trans_results.txt" \
    --trans-p-value "$TRANS_PVAL" \
    # --no-fdr

echo "  >> Cis Results:   $OUTPUT_DIR/cis_results.txt"
echo "  >> Trans Results: $OUTPUT_DIR/trans_results.txt"


# --- Mode B: Location-Free Analysis (All-Pairs) ---
# Uncomment the block below to run this.
# This ignores location files and tests ALL gene-SNP pairs.
# Useful if you don't have location data or want a pure ANOVA-style analysis.
# Warning: Output can be huge.

# echo "-----------------------------------------------------------------------------"
# echo "[Step 5b] Running All-Pairs Analysis (Location-Free)..."
#
# python3 "$CODE_DIR/src/matrix_eqtl/analysis/run_matrix_eqtl.py" \
#     --genotype-matrix "$OUTPUT_DIR/preprocessing/genotype.matrix" \
#     --gene-expression-matrix "$OUTPUT_DIR/preprocessing/expression.qnorm" \
#     --covariates "$OUTPUT_DIR/covariates/final_covariates.txt" \
#     --output-file "$OUTPUT_DIR/all_pairs_results.txt" \
#     --p-value "$ALL_PAIRS_PVAL" \
#     --no-fdr
#
# echo "  >> All-Pairs Results: $OUTPUT_DIR/all_pairs_results.txt"

echo "-----------------------------------------------------------------------------"
echo ">>> PIPELINE COMPLETED SUCCESSFULLY!"
echo "-----------------------------------------------------------------------------"
