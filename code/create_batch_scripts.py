import pandas as pd
import argparse
import os
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Generate MatrixEQTL scripts from CSV")
    parser.add_argument('--csv', required=True, help='Input CSV file with parameters')
    parser.add_argument('--code-dir', required=True, help='Path to the code directory containing scripts')
    return parser.parse_args()

def generate_matrix_script(row, code_dir):
    job_id = str(row['id'])
    work_dir = str(row['work_dir'])
    peer_n = int(row.get('peer', 15))
    cis_dist = int(row.get('cis_dist', 1000000))
    p_cis = float(row.get('p_cis', 1e-5))
    p_trans = float(row.get('p_trans', 1e-5))
    maf = float(row.get('maf', 0.05))
    threads = int(row.get('threads', 4))
    vcf_file = str(row['vcf'])
    exp_file = str(row['expression'])
    
    has_loc = False
    loc_file = ""
    if 'location' in row and pd.notna(row['location']):
        loc_str = str(row['location']).strip()
        if loc_str and loc_str.lower() != 'nan' and loc_str.lower() != 'na':
            has_loc = True
            loc_file = loc_str

    script = f"""#!/bin/bash
# Job ID: {job_id} (MatrixEQTL)
export WORK_DIR="{work_dir}"
export CODE_DIR="{code_dir}"
# Separate output dir for MatrixEQTL
export OUT_DIR="${{WORK_DIR}}/output_matrix" 

export VCF_FILE="{vcf_file}"
export EXP_FILE="{exp_file}"
export THREADS={threads}

mkdir -p ${{OUT_DIR}}/step1_plink
mkdir -p ${{OUT_DIR}}/step2_covariates
mkdir -p ${{OUT_DIR}}/final_cis
# Logs
exec > >(tee -i ${{WORK_DIR}}/run_matrix.log) 2>&1
echo ">>> Starting MatrixEQTL Job {job_id}"

# Step 1: VCF (PLINK)
echo ">>> Step 1: Processing VCF for MatrixEQTL..."
python3 ${{CODE_DIR}}/process_vcf_plink.py --vcf "${{VCF_FILE}}" --expression "${{EXP_FILE}}" --out-dir "${{OUT_DIR}}/step1_plink" --maf {maf} --threads ${{THREADS}} --snps-only

# Step 2: Covariates
echo ">>> Step 2: Generating Covariates..."
python3 ${{CODE_DIR}}/run_iqn.py -i "${{EXP_FILE}}" -o "${{OUT_DIR}}/step2_covariates/expression.qnorm"
# PCA
plink --bfile "${{OUT_DIR}}/step1_plink/plink_temp" --pca 3 --threads ${{THREADS}} --out "${{OUT_DIR}}/step2_covariates/plink_pca"
python3 ${{CODE_DIR}}/format_plink_pca.py -i "${{OUT_DIR}}/step2_covariates/plink_pca.eigenvec" -o "${{OUT_DIR}}/step2_covariates/genotype.pcs" -n 3
# PEER
python3 ${{CODE_DIR}}/run_peer.py -i "${{OUT_DIR}}/step2_covariates/expression.qnorm" -n {peer_n} -o "${{OUT_DIR}}/step2_covariates/peer_factors.tsv"
# Combine
python3 ${{CODE_DIR}}/combine_covariates.py -p "${{OUT_DIR}}/step2_covariates/genotype.pcs" -f "${{OUT_DIR}}/step2_covariates/peer_factors.tsv" -o "${{OUT_DIR}}/step2_covariates/final_covariates.txt"

# Step 3: Analysis
"""
    if has_loc:
        script += f"""
echo ">>> Step 3: Running MatrixEQTL (Standard)..."
export GENE_LOC_FILE="${{OUT_DIR}}/gene_locations_fixed.tsv"
if [ ! -f "${{GENE_LOC_FILE}}" ]; then sed '1s/features/geneid/' "{loc_file}" > "${{GENE_LOC_FILE}}"; fi
python3 ${{CODE_DIR}}/run_matrix_eqtl.py --genotype-matrix "${{OUT_DIR}}/step1_plink/genotype.matrix" --genotype-positions "${{OUT_DIR}}/step1_plink/snp_positions.txt" --gene-expression-matrix "${{OUT_DIR}}/step2_covariates/expression.qnorm" --gene-positions "${{GENE_LOC_FILE}}" --covariates "${{OUT_DIR}}/step2_covariates/final_covariates.txt" --output-file "${{OUT_DIR}}/final_cis/cis_results.txt" --trans-output-file "${{OUT_DIR}}/final_cis/trans_results.txt" --p-value {p_cis} --trans-p-value {p_trans} --cis-distance {cis_dist} --chunk-size 2000
echo "Done. Results in ${{OUT_DIR}}/final_cis/"
"""
    else:
        script += f"""
echo ">>> Step 3: Running MatrixEQTL (Location-Free)..."
python3 ${{CODE_DIR}}/run_matrix_eqtl_noloc.py --genotype-matrix "${{OUT_DIR}}/step1_plink/genotype.matrix" --expression-matrix "${{OUT_DIR}}/step2_covariates/expression.qnorm" --covariates "${{OUT_DIR}}/step2_covariates/final_covariates.txt" --output-file "${{OUT_DIR}}/final_cis/global_results.txt" --p-value {p_trans} --no-fdr
echo "Done. Results in ${{OUT_DIR}}/final_cis/"
"""
    return script

def generate_qtltools_script(row, code_dir):
    job_id = str(row['id'])
    work_dir = str(row['work_dir'])
    peer_n = int(row.get('peer', 15))
    cis_dist = int(row.get('cis_dist', 1000000))
    maf = float(row.get('maf', 0.05))
    threads = int(row.get('threads', 4))
    vcf_file = str(row['vcf'])
    exp_file = str(row['expression'])
    
    # ID Type (new optional column)
    id_type = str(row.get('id_type', 'ensembl_gene_id')).strip()
    
    # Check location logic
    # Options:
    # 1. Path to file -> Use --positions
    # 2. "fetch" or "auto" -> Use --fetch-pos
    # 3. Empty/NaN -> Skip
    
    use_fetch = False
    loc_file = ""
    
    if 'location' in row and pd.notna(row['location']):
        val = str(row['location']).strip()
        if val.lower() in ['fetch', 'auto', 'biomart']:
            use_fetch = True
        elif val and val.lower() != 'nan':
            loc_file = val
    
    # If neither file nor fetch, cannot run QTLtools cis
    if not loc_file and not use_fetch:
        return None

    script = f"""#!/bin/bash
# Job ID: {job_id} (QTLtools)
export WORK_DIR="{work_dir}"
export CODE_DIR="{code_dir}"
# Separate output dir for QTLtools
export OUT_DIR="${{WORK_DIR}}/output_qtltools"

export VCF_FILE="{vcf_file}"
export EXP_FILE="{exp_file}"
export THREADS={threads}

mkdir -p ${{OUT_DIR}}/step1_plink
mkdir -p ${{OUT_DIR}}/step2_covariates
mkdir -p ${{OUT_DIR}}/qtltools_input
mkdir -p ${{OUT_DIR}}/qtltools_results
exec > >(tee -i ${{WORK_DIR}}/run_qtltools.log) 2>&1
echo ">>> Starting QTLtools Job {job_id}"

# Note: We need to run Steps 1 & 2 to generate Covariates and ensure sample matching
# Step 1: VCF (PLINK) - Needed for PCA and sample list
echo ">>> Step 1: Processing VCF (for PCA/Covariates)..."
python3 ${{CODE_DIR}}/process_vcf_plink.py --vcf "${{VCF_FILE}}" --expression "${{EXP_FILE}}" --out-dir "${{OUT_DIR}}/step1_plink" --maf {maf} --threads ${{THREADS}} --snps-only

# Step 2: Covariates
echo ">>> Step 2: Generating Covariates..."
python3 ${{CODE_DIR}}/run_iqn.py -i "${{EXP_FILE}}" -o "${{OUT_DIR}}/step2_covariates/expression.qnorm"

# PCA
plink --bfile "${{OUT_DIR}}/step1_plink/plink_temp" --pca 3 --threads ${{THREADS}} --out "${{OUT_DIR}}/step2_covariates/plink_pca"
python3 ${{CODE_DIR}}/format_plink_pca.py -i "${{OUT_DIR}}/step2_covariates/plink_pca.eigenvec" -o "${{OUT_DIR}}/step2_covariates/genotype.pcs" -n 3
# PEER
python3 ${{CODE_DIR}}/run_peer.py -i "${{OUT_DIR}}/step2_covariates/expression.qnorm" -n {peer_n} -o "${{OUT_DIR}}/step2_covariates/peer_factors.tsv"
# Combine
python3 ${{CODE_DIR}}/combine_covariates.py -p "${{OUT_DIR}}/step2_covariates/genotype.pcs" -f "${{OUT_DIR}}/step2_covariates/peer_factors.tsv" -o "${{OUT_DIR}}/step2_covariates/final_covariates.txt"


echo ">>> Step 3: Preparing QTLtools Input..."
"""

    if use_fetch:
        script += f"""
# Mode: Auto-Fetch from BioMart
echo "    Fetching positions for ID type: {id_type}..."
python3 ${{CODE_DIR}}/qtltools_prep.py \\
    --expression "${{OUT_DIR}}/step2_covariates/expression.qnorm" \\
    --fetch-pos \\
    --id-type "{id_type}" \\
    --out "${{OUT_DIR}}/qtltools_input/expression.bed"
"""
    else:
        script += f"""
# Mode: Local Position File
export GENE_LOC_FILE="${{OUT_DIR}}/gene_locations_fixed.tsv"
if [ ! -f "${{GENE_LOC_FILE}}" ]; then sed '1s/features/geneid/' "{loc_file}" > "${{GENE_LOC_FILE}}"; fi

python3 ${{CODE_DIR}}/qtltools_prep.py \\
    --expression "${{OUT_DIR}}/step2_covariates/expression.qnorm" \\
    --positions "${{GENE_LOC_FILE}}" \\
    --out "${{OUT_DIR}}/qtltools_input/expression.bed"
"""

    script += f"""
# Prepare Covariates
cp "${{OUT_DIR}}/step2_covariates/final_covariates.txt" "${{OUT_DIR}}/qtltools_input/covariates.txt"

echo ">>> Step 4: Running QTLtools (Permutation)..."
QTLtools cis \\
    --vcf "${{VCF_FILE}}" \\
    --bed "${{OUT_DIR}}/qtltools_input/expression.bed.gz" \\
    --cov "${{OUT_DIR}}/qtltools_input/covariates.txt" \\
    --out "${{OUT_DIR}}/qtltools_results/permutations.txt" \\
    --permute 1000 \\
    --window {cis_dist} \\
    --normal

echo "Done. Results in ${{OUT_DIR}}/qtltools_results/permutations.txt"
"""
    return script

def main():
    args = get_args()
    try:
        df = pd.read_csv(args.csv)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)
        
    scripts_dir = "generated_scripts"
    if not os.path.exists(scripts_dir):
        os.makedirs(scripts_dir)
        
    master_script = open("run_all_jobs.sh", "w")
    master_script.write("#!/bin/bash\n")
    
    for idx, row in df.iterrows():
        job_id = str(row['id'])
        
        # 1. Matrix Script
        m_content = generate_matrix_script(row, os.path.abspath(args.code_dir))
        m_name = f"run_{job_id}_matrix.sh"
        with open(os.path.join(scripts_dir, m_name), 'w') as f:
            f.write(m_content)
            
        # 2. QTLtools Script (Only if location exists)
        q_content = generate_qtltools_script(row, os.path.abspath(args.code_dir))
        q_name = f"run_{job_id}_qtltools.sh"
        if q_content:
            with open(os.path.join(scripts_dir, q_name), 'w') as f:
                f.write(q_content)
            print(f"Generated: {m_name} & {q_name}")
            master_script.write(f"bash {scripts_dir}/{m_name}\n")
            master_script.write(f"bash {scripts_dir}/{q_name}\n")
        else:
            print(f"Generated: {m_name} (No QTLtools due to missing location)")
            master_script.write(f"bash {scripts_dir}/{m_name}\n")
            
    master_script.close()
    print(f"Scripts saved to {scripts_dir}/")


if __name__ == '__main__':
    main()
