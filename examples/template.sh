#!/bin/bash
# Wrapper to run the batch script generator.
# Usage: bash template.sh <csv_file>

if [ -z "$1" ]; then
    echo "Usage: bash template.sh <csv_file>"
    exit 1
fi

CSV_FILE=$1
CODE_DIR="$(dirname $0)/matrix_eqtl_pipeline/code" # Assuming script is in root

# If CODE_DIR doesn't match the current structure, try to find it
if [ ! -d "$CODE_DIR" ]; then
    # Maybe we are inside matrix_eqtl_pipeline root?
    if [ -d "code" ]; then
        CODE_DIR="$(pwd)/code"
    elif [ -d "/data/work/new_QTL/matrix_eqtl_pipeline/code" ]; then
       # Hardcoded fallback based on user context
       CODE_DIR="/data/work/new_QTL/matrix_eqtl_pipeline/code"
    else
       echo "Error: Could not locate 'code' directory."
       exit 1
    fi
fi

echo "Using Code Dir: $CODE_DIR"
python3 ${CODE_DIR}/create_batch_scripts.py --csv "$CSV_FILE" --code-dir "$CODE_DIR"

chmod +x generated_scripts/*.sh
chmod +x run_all_jobs.sh

echo ""
echo ">>> Generation Complete."
echo "You can now run individual scripts in 'generated_scripts/'"
echo "Or run 'bash run_all_jobs.sh' to execute them all sequentially."
