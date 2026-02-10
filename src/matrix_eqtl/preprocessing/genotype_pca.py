import argparse
import subprocess
from pathlib import Path
import pandas as pd
import os
import sys

def run_pca(genotype_prefix, output_dir, pca_n=3, threads=1, plink_bin="plink2", plink_version=2):
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
        plink_bin, 
        "--bfile", str(genotype_prefix), 
        "--pca", str(pca_n), 
        "--out", str(pca_out),
        "--threads", str(threads)
    ]

    if plink_version == 2:
        # PLINK 2 optimization: approx defaults to fast
        # Ensure it works same as v1.9 basic PCA
        cmd.append("--bad-freqs") # Allow small sample sizes (<50) without error
    else:
        # PLINK 1.9
        pass
    
    print(f"Running PLINK PCA ({plink_bin} v{plink_version}): {' '.join(cmd)}")
    subprocess.check_call(cmd)
    
    # Format for MatrixEQTL (Transpose)
    vec_file = str(pca_out) + ".eigenvec"
    if not os.path.exists(vec_file):
        print("ERROR: PLINK did not produce .eigenvec file.")
        sys.exit(1)
        
    df_pca = pd.read_csv(vec_file, sep=r'\s+', header=None)
    
    if plink_version == 2:
        # PLINK 2 output format: FID IID PC1... (header is usually commented #FID IID)
        # But pandas read_csv with header=None reads it as data if comments not handled
        # Actually plink2 output often has a header line starting with #
        # Let's re-read with comment handling
        df_pca = pd.read_csv(vec_file, sep=r'\s+', comment='#', header=None)
        # Columns: [FID aka IID], [IID], PC1... ?
        # PLINK 2 defaults: FID IID PC1...
        # Since we use --double-id, FID=IID.
        # So col 0 = IID, col 1 = IID, col 2 = PC1...
        # We need IID and PCs
        df_pca = df_pca.iloc[:, 1:] 
    else:
        # PLINK 1.9: FID IID PC1...
        df_pca = df_pca.iloc[:, 1:] 

    # Rename columns: IID, PC1, ...
    cols = ['id'] + [f'PC{i+1}' for i in range(pca_n)]
    # Safety check on column count
    if len(df_pca.columns) != len(cols):
        # Maybe PLINK outputted more or fewer cols?
        # Just slice what we need
        df_pca = df_pca.iloc[:, :len(cols)]
        
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
    
    # PLINK Version Control
    parser.add_argument("--plink-bin", default="plink2", help="Path/name of PLINK executable")
    parser.add_argument("--plink-version", type=int, default=2, choices=[1, 2], help="PLINK version (1 or 2)")
    
    args = parser.parse_args()
    
    print(f"DEBUG: Python CWD: {os.getcwd()}")
    
    run_pca(args.genotype, args.out_dir, args.pca_n, args.threads, 
            plink_bin=args.plink_bin, plink_version=args.plink_version)

if __name__ == "__main__":
    main()
