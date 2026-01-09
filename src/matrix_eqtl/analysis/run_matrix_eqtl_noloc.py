from __future__ import print_function
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import sys
import argparse
import os
from check_file import check_file

def get_args():
    parser = argparse.ArgumentParser(description="Run MatrixEQTL in location-agnostic mode (Engine)")
    parser.add_argument('-m', '--genotype-matrix', required=True, help='Genotype matrix file')
    parser.add_argument('-e', '--expression-matrix', required=True, help='Phenotype (Expression/Metabolite) matrix file')
    parser.add_argument('-c', '--covariates', help='Covariates file (optional)')
    parser.add_argument('-o', '--output-file', default='MatrixEqtlOutput_NoLoc.txt', help='Output file path')
    parser.add_argument('-v', '--p-value', type=float, default=1e-5, help='P-value threshold')
    parser.add_argument('--no-fdr', action='store_true', help='Do not calculate FDR (saves memory)')
    parser.add_argument('--chunk-size', type=int, default=2000, help='Slice size')
    
    args = parser.parse_args()

    check_file(args.genotype_matrix)
    check_file(args.expression_matrix)
    if args.covariates:
        check_file(args.covariates)

    return args

def main():
    args = get_args()

    r = ro.r
    
    # Source the new R script
    # Assumes mxeqtl_noloc.R is in the code/ directory relative to execution or fixed path
    # We will try to find it relative to this script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    r_script_path = os.path.join(script_dir, "mxeqtl_noloc.R")
    
    if not os.path.exists(r_script_path):
        print(f"Error: Could not find R script at {r_script_path}")
        sys.exit(1)

    r.source(r_script_path)

    # Call the run_engine function
    # Arguments: snp_file, expr_file, cov_file, output_file, p_value, no_fdr, chunk_size
    print(f"Running Matrix_eQTL_engine...")
    print(f"  Genotypes: {args.genotype_matrix}")
    print(f"  Phenotypes: {args.expression_matrix}")
    print(f"  Covariates: {args.covariates}")
    print(f"  Output: {args.output_file}")
    print(f"  P-value threshold: {args.p_value}")
    print(f"  No FDR: {args.no_fdr}")
    print(f"  Chunk Size: {args.chunk_size}")
    
    cov_arg = args.covariates if args.covariates else ""
    
    r.run_engine(
        args.genotype_matrix, 
        args.expression_matrix, 
        cov_arg, 
        args.output_file, 
        args.p_value,
        args.no_fdr,
        args.chunk_size
    )
    
    print("Done!")

if __name__ == '__main__':
    main()
