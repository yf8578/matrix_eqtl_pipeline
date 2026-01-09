#!/bin/bash

# =============================================================================
# MatrixEQTL Pipeline - Custom Analysis Script (Robust Version)
# =============================================================================

# Add error handling: Stop on first error
set -e 

# 1. Set paths
PROJ_DIR="/data/work/QTLpipeline/matrix_eqtl_pipeline"
cd $PROJ_DIR
export PYTHONPATH=$PROJ_DIR/src

VCF_FILE="/data/work/SNP_data/control/only_control_20240613_v2.vcf.gz"
EXPR_FILE="/data/work/xQTL/pQTL/01_data/control/49_control_protein_expression_matrix.tsv"
OUTPUT_DIR="results_control"

mkdir -p $OUTPUT_DIR



# -----------------------------------------------------------------------------
# Step 0: Gene Locations (Essential for Cis-analysis)
# -----------------------------------------------------------------------------
echo "[Step 0] Generating Gene Locations..."
# Note: Added check to not fail if network is down, just warn.
set +e
python3 src/matrix_eqtl/preprocessing/locations.py \
    --expression "$EXPR_FILE" \
    --out "$OUTPUT_DIR/gene_locations.txt" \
    --format matrixeqtl \
    --fetch-pos \
    --id-type ensembl_gene_id
    
if [ $? -ne 0 ]; then
    echo "Warning: Step 0 failed (probably network). Continuing..."
else
    echo "Step 0 Done."
fi
set -e

# -----------------------------------------------------------------------------
# Step 1: Intersection
# -----------------------------------------------------------------------------
echo "[Step 1] Intersecting Samples..."
python3 src/matrix_eqtl/utils/io.py \
    --vcf "$VCF_FILE" \
    --expression "$EXPR_FILE" \
    --out "$OUTPUT_DIR/samples.txt"

# -----------------------------------------------------------------------------
# Step 2: Expression (Protein)
# -----------------------------------------------------------------------------
echo "[Step 2] Processing Expression..."
python3 src/matrix_eqtl/preprocessing/expression.py \
    --expression "$EXPR_FILE" \
    --samples "$OUTPUT_DIR/samples.txt" \
    --out-dir "$OUTPUT_DIR/preprocessing" \
    --no-tpm-filter \
    --max-missing 0.1

# -----------------------------------------------------------------------------
# Step 3: Genotypes
# -----------------------------------------------------------------------------
echo "[Step 3] Processing Genotypes..."
python3 src/matrix_eqtl/preprocessing/genotype.py \
    --vcf "$VCF_FILE" \
    --samples "$OUTPUT_DIR/samples.txt" \
    --out-dir "$OUTPUT_DIR/preprocessing" \
    --tool matrixeqtl \
    --maf 0.01 --geno 0.05 --hwe 1e-6 \
    --threads 5

# -----------------------------------------------------------------------------
# Step 4: Covariates
# -----------------------------------------------------------------------------
echo "[Step 4] Covariates..."
# Step 4a: Genotype PCA
# -----------------------------------------------------------------------------
echo "[Step 4a] Genotype PCA..."
python3 src/matrix_eqtl/preprocessing/genotype_pca.py \
    --genotype "$OUTPUT_DIR/preprocessing/genotypes" \
    --out-dir "$OUTPUT_DIR/covariates" \
    --pca-n 3 \
    --threads 4

# -----------------------------------------------------------------------------
# Step 4b: PEER Factors
# -----------------------------------------------------------------------------
echo "[Step 4b] PEER Factors..."
python3 src/matrix_eqtl/preprocessing/expression_peer.py \
    --expression "$OUTPUT_DIR/preprocessing/expression.qnorm" \
    --out-dir "$OUTPUT_DIR/covariates" \
    --peer-n 5

# -----------------------------------------------------------------------------
# Step 4c: Combine Covariates
# -----------------------------------------------------------------------------
echo "[Step 4c] Combining Covariates..."
# Concatenate the files (Logic: PCA + PEER)
# We can use a simple python oneliner or cat if we handle headers carefully.
# Let's simple cat them, assuming both have headers and same columns?
# Actually headers might conflict.
# Better to use a small merger script or reuse covariates.py just for merging if arguments allow?
# `covariates.py` does logic.
# Let's just create a simple merge script or append.
# Or just use the original covariates.py but with pca=0 and peer=0? 
# No, `covariates.py` generates them.

# Let's perform a robust merge using python in-line
python3 -c "
import pandas as pd
import os
out_dir = '$OUTPUT_DIR/covariates'
pca_file = os.path.join(out_dir, 'genotype_pcs.txt')
peer_file = os.path.join(out_dir, 'peer_factors.tsv')
out_file = os.path.join(out_dir, 'final_covariates.txt')

dfs = []
if os.path.exists(pca_file):
    dfs.append(pd.read_csv(pca_file, sep='\t', index_col=0))
if os.path.exists(peer_file):
    dfs.append(pd.read_csv(peer_file, sep='\t', index_col=0))

if dfs:
    final = pd.concat(dfs, axis=0, join='inner')
    final.to_csv(out_file, sep='\t')
    print(f'Combined covariates saved to {out_file}')
else:
    print('No covariates found to combine.')
" \
    --threads 4

# -----------------------------------------------------------------------------
# Step 5: Run
# -----------------------------------------------------------------------------
echo "[Step 5] Analysis..."
# Check if location file exists (from step 0)
POS_ARGS=""
if [ -f "$OUTPUT_DIR/gene_locations.txt" ]; then
    POS_ARGS="--gene-positions $OUTPUT_DIR/gene_locations.txt --genotype-positions $OUTPUT_DIR/preprocessing/snp_positions.txt"
fi

python3 src/matrix_eqtl/analysis/run_matrix_eqtl.py \
    --genotype-matrix "$OUTPUT_DIR/preprocessing/genotype.matrix" \
    --gene-expression-matrix "$OUTPUT_DIR/preprocessing/expression.qnorm" \
    --covariates "$OUTPUT_DIR/covariates/final_covariates.txt" \
    $POS_ARGS \
    --output-file "$OUTPUT_DIR/cis_results.txt" \
    --cis-distance 1000000 
    
echo ">>> SUCCESS! Pipeline finished."
