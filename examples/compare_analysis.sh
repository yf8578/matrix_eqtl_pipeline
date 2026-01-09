#!/bin/bash
# compare_analysis.sh
# Run MatrixEQTL in 3 Modes to Compare Results

export WORK_DIR=${WORK_DIR:-/data/work/new_QTL}
export CODE_DIR=${WORK_DIR}/matrix_eqtl_pipeline/code
# Or default if user didn't set env:
# export CODE_DIR=${CODE_DIR:-/data/work/new_QTL/code}
export DEMO_DIR=${WORK_DIR}/demo_data
export OUT_DIR=${DEMO_DIR}/comparison_results

mkdir -p ${OUT_DIR}

COMMON_ARGS="--chunk-size 1000 --sep \t"

echo "=========================================================="
echo " Comparison 1: Standard Mode (Distinguish Cis vs Trans)"
echo " Cis-Distance = 1Mb. P-value < 0.05"
echo "=========================================================="

python3 ${CODE_DIR}/run_matrix_eqtl.py \
    --genotype-matrix ${DEMO_DIR}/genotype.demo.txt \
    --genotype-positions ${DEMO_DIR}/snp_pos.demo.txt \
    --gene-expression-matrix ${DEMO_DIR}/expression.demo.txt \
    --gene-positions ${DEMO_DIR}/gene_pos.demo.txt \
    --covariates ${DEMO_DIR}/covariates.demo.txt \
    --output-file ${OUT_DIR}/standard_cis.txt \
    --trans-output-file ${OUT_DIR}/standard_trans.txt \
    --p-value 0.05 \
    --trans-p-value 0.05 \
    --cis-distance 1000000 \
    ${COMMON_ARGS}

echo "  -> Cis Results: ${OUT_DIR}/standard_cis.txt"
echo "  -> Trans Results: ${OUT_DIR}/standard_trans.txt"


echo "=========================================================="
echo " Comparison 2: Classic Mode (Location-Free / All-vs-All)"
echo " No Position Files Provided. Calculating ALL pairs."
echo " P-value < 0.05"
echo "=========================================================="

python3 ${CODE_DIR}/run_matrix_eqtl.py \
    --genotype-matrix ${DEMO_DIR}/genotype.demo.txt \
    --gene-expression-matrix ${DEMO_DIR}/expression.demo.txt \
    --covariates ${DEMO_DIR}/covariates.demo.txt \
    --output-file /dev/null \
    --trans-output-file ${OUT_DIR}/classic_all.txt \
    --trans-p-value 0.05 \
    --p-value 0 \
    --cis-distance 0 \
    ${COMMON_ARGS}

echo "  -> Classic Results: ${OUT_DIR}/classic_all.txt"

echo "=========================================================="
echo " Summary of Counts (P < 0.05)"
echo "=========================================================="
echo "Standard Cis:"
wc -l ${OUT_DIR}/standard_cis.txt
echo "Standard Trans:"
wc -l ${OUT_DIR}/standard_trans.txt
echo "Classic (All-vs-All):"
wc -l ${OUT_DIR}/classic_all.txt

echo "Done! Check ${OUT_DIR} for files."
