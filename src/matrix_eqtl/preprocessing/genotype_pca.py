import argparse
import subprocess
from pathlib import Path
import pandas as pd
import os
import sys

def run_pca(genotype_prefix, output_dir, pca_n=3, threads=1):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check input
    bed_file = Path(str(genotype_prefix) + ".bed")
    if not bed_file.exists():
        print(f"ERROR: Genotype BED file not found at: {bed_file.resolve()}")
        # Debug helper
        parent_dir = bed_file.parent
        if parent_dir.exists():
            print(f"Contents of {parent_dir}:")
            for f in parent_dir.iterdir():
                print(f"  {f.name}")
        else:
             print(f"Directory {parent_dir} NOT found.")
        sys.exit(1)

    pca_out = output_dir / "genotype_pca"
    
    # Run PLINK
    cmd = [
        "plink", 
        "--bfile", str(genotype_prefix), 
        "--pca", str(pca_n), 
        "--out", str(pca_out),
        "--threads", str(threads)
    ]
    
    print(f"Running PLINK PCA: {' '.join(cmd)}")
    subprocess.check_call(cmd)
    
    # Format for MatrixEQTL (Transpose)
    vec_file = str(pca_out) + ".eigenvec"
    if not os.path.exists(vec_file):
        print("ERROR: PLINK did not produce .eigenvec file.")
        sys.exit(1)
        
    df_pca = pd.read_csv(vec_file, sep=r'\s+', header=None)
    # FID IID PC1 PC2...
    # Keep IID and PCs
    # Drop FID (col 0)
    df_pca = df_pca.iloc[:, 1:] 
    # Rename columns: IID, PC1, ...
    cols = ['id'] + [f'PC{i+1}' for i in range(pca_n)]
    df_pca.columns = cols
    
    # Transpose so rows=PCs, cols=Samples
    df_pca = df_pca.set_index('id').T
    
    final_out = output_dir / "genotype_pcs.txt"
    df_pca.to_csv(final_out, sep='\t')
    print(f"PCA Covariates saved to: {final_out}")

def main():
    parser = argparse.ArgumentParser(description="Genotype PCA (PLINK)")
    parser.add_argument("--genotype", required=True, help="PLINK prefix")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--pca-n", type=int, default=3, help="Number of PCs")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    
    args = parser.parse_args()
    
    print(f"DEBUG: Python CWD: {os.getcwd()}")
    
    run_pca(args.genotype, args.out_dir, args.pca_n, args.threads)

if __name__ == "__main__":
    main()
