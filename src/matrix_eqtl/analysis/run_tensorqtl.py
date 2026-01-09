
import sys
import os
import argparse
import pandas as pd
import logging

# Add local tensorqtl to path
current_dir = os.path.dirname(os.path.abspath(__file__))
tensorqtl_dir = os.path.join(current_dir, 'tensorqtl-1.0.10')
sys.path.insert(0, tensorqtl_dir)

import tensorqtl
from tensorqtl import genotypeio, cis, trans

def get_args():
    parser = argparse.ArgumentParser(description="Run TensorQTL wrapper")
    parser.add_argument('--plink', required=True, help='PLINK genotype prefix')
    parser.add_argument('--phenotypes', required=True, help='Phenotype BED file')
    parser.add_argument('--covariates', required=True, help='Covariates text file')
    parser.add_argument('--out', required=True, help='Output prefix')
    parser.add_argument('--mode', default='cis', choices=['cis', 'cis_nominal', 'trans'], help='Analysis mode')
    parser.add_argument('--window', type=int, default=1000000, help='Window size (default 1Mb)')
    parser.add_argument('--fdr', type=float, default=0.05, help='FDR threshold for cis-QTL perm')
    return parser.parse_args()

def main():
    args = get_args()
    
    # Check GPU availability (just logging)
    import torch
    if torch.cuda.is_available():
        print(f"TensorQTL: Using GPU {torch.cuda.get_device_name(0)}")
    else:
        print("TensorQTL: Running on CPU (Warning: slower)")

    # Load Phenotypes
    print(f"Loading phenotypes: {args.phenotypes}")
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.phenotypes)
    
    # Load Covariates
    print(f"Loading covariates: {args.covariates}")
    # Covariates file is usually Covariate x Sample for QTLtools/MatrixEQTL
    # TensorQTL expects Sample x Covariate usually when reading from file, OR dataframe
    # The README says: "Tab-delimited text file (covariates x samples) or dataframe (samples x covariates)"
    # BUT standard Pandas read is Sample x Feature usually.
    # We transpose to match expected dataframe structure (samples x covariates) if input is rows=covariates
    cov = pd.read_csv(args.covariates, sep='\t', index_col=0)
    covariates_df = cov.T 
    
    # Load Genotypes
    print(f"Loading genotypes: {args.plink}")
    pr = genotypeio.PlinkReader(args.plink)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    
    # Intersect Samples
    # tensorqtl handles this but let's be explicit
    common_samples = sorted(list(set(covariates_df.index) & set(phenotype_df.columns) & set(genotype_df.columns)))
    print(f"Matched samples: {len(common_samples)}")
    
    if len(common_samples) == 0:
        print("Error: No overlapping samples found.")
        sys.exit(1)
        
    # Run Mode
    if args.mode == 'cis':
        print(f"Running Cis-QTL Permutations (Window={args.window})...")
        cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, 
                             covariates_df, window=args.window, seed=123456)
        
        # Calculate Q-values
        tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
        
        out_file = f"{args.out}.cis_perm.txt.gz"
        cis_df.to_csv(out_file, sep='\t', compression='gzip')
        print(f"Saved permutation results: {out_file}")
        
    elif args.mode == 'cis_nominal':
        print(f"Running Cis-QTL Nominal (Window={args.window})...")
        # output_dir is derived from args.out prefix
        out_dir = os.path.dirname(args.out)
        prefix_base = os.path.basename(args.out)
        
        cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                        prefix_base, covariates_df, window=args.window, 
                        output_dir=out_dir)
        print(f"Saved nominal parquet files to {out_dir}")
        
    elif args.mode == 'trans':
        print("Running Trans-QTL...")
        trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df,
                                   return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
        # Filter cis
        # trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)
        out_file = f"{args.out}.trans.parquet"
        trans_df.to_parquet(out_file)
        print(f"Saved trans results: {out_file}")

if __name__ == "__main__":
    main()
