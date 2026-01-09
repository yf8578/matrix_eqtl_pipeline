from __future__ import print_function
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import sys
import io
import argparse
import os
from check_file import check_file


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--genotype-matrix', required=True, help='')
    parser.add_argument('-p', '--genotype-positions', required=False, default=None, help='')
    parser.add_argument('-e', '--gene-expression-matrix', required=True, help='')
    parser.add_argument('-g', '--gene-positions', required=False, default=None, help='')
    parser.add_argument('-c', '--covariates', help='')
    parser.add_argument('-o', '--output-file', default='MatrixEqtlOutput', help='')
    parser.add_argument('-v', '--p-value', type=float, default=0.05, help='')
    parser.add_argument('-q', '--qq-plot', default='MatrixEqtlQQPlot.pdf', help='')
    parser.add_argument('--trans-output-file', default="", help='')
    parser.add_argument('--trans-p-value', type=float, default=0.0, help='')
    parser.add_argument('--model', default='linear', choices={'linear', 'anova', 'linear_cross'}, help='')
    parser.add_argument('--cis-distance', type=float, default=1e6, help='')
    parser.add_argument('--maf', default=0.0, help='')
    parser.add_argument('--no-header', action='store_false', help='')
    parser.add_argument('--no-rownames', action='store_false', help='')
    parser.add_argument('--missing', default='NA', help='')
    parser.add_argument('--sep', default='\t', help='')
    parser.add_argument('--chunk-size', type=int, default=2000, help='')
    parser.add_argument('--no-fdr', action='store_true', help='Do not calculate FDR (saves memory)')

    args = parser.parse_args()

    check_file(args.genotype_matrix)

    return(args)


import pandas as pd

def align_inputs(args):
    """
    Check if sample columns in Expression and Covariates match Genotype.
    If not, create temporary aligned files and update args to point to them.
    Genotype is considered the 'Anchor/Master'.
    """
    print(">>> Checking Sample Alignment...")
    
    # helper to get samples
    def get_samples(fpath, sep):
        with open(fpath, 'r') as f:
            header = f.readline().strip().split(sep)
            # Assuming first col is ID
            return header[1:]

    sep_map = {'\\t': '\t', ',': ','} # handle escaped tab if passed
    sep = sep_map.get(args.sep, args.sep)

    # 1. Get Master Samples from Genotype
    try:
        master_samples = get_samples(args.genotype_matrix, sep)
    except Exception as e:
        print(f"Error reading genotype header: {e}")
        return # Let downstream fail if file is bad

    # 2. Check Expression
    expr_samples = get_samples(args.gene_expression_matrix, sep)
    
    # 3. Check Covariates (if present)
    cov_samples = None
    if args.covariates:
        cov_samples = get_samples(args.covariates, sep)

    # Validate Intersection
    # Should we intersect? Or assume Genotype is correct?
    # Usually Genotype is from PLINK, so samples might be dropped.
    # The safest is: Align to Genotype. If Expr/Cov missing Genotype samples -> Error.
    
    need_align_expr = (expr_samples != master_samples)
    need_align_cov = (cov_samples is not None and cov_samples != master_samples)

    if not need_align_expr and not need_align_cov:
        print("  ✅ All files are already aligned.")
        return

    print("  ⚠️ Sample mismatch detected. Aligning files to Genotype order (In-Place)...")
    
    # Function to rewrite file in-place
    def rewrite_in_place(in_path, name):
        print(f"    Aligning {name} (Overwriting original)...")
        # Optimization: Read only header first to confirm? No, we need to read all to write.
        # But we do it carefully.
        
        try:
             df = pd.read_csv(in_path, sep=sep)
        except Exception as e:
             print(f"Error reading {name}: {e}")
             sys.exit(1)
             
        id_col = df.columns[0]
        
        # Check missing (Genotype samples must exist in Target)
        cols = df.columns.tolist()
        missing = [s for s in master_samples if s not in cols]
        if missing:
            print(f"    Error: Genotype samples {missing[:5]}... not found in {name} file.")
            sys.exit(1)
            
        # Reorder
        new_order = [id_col] + list(master_samples)
        df_aligned = df[new_order]
        
        # Write to temp then move
        tmp_path = in_path + ".tmp"
        df_aligned.to_csv(tmp_path, sep=sep, index=False, float_format='%.5g')
        os.replace(tmp_path, in_path)
        print(f"    ✅ Updated {in_path}")

    # Perform Alignment
    if need_align_expr:
        rewrite_in_place(args.gene_expression_matrix, "Expression")
        
    if need_align_cov:
        rewrite_in_place(args.covariates, "Covariates")
        
    print("  ✅ Alignment passed. Files are now consistent with Genotype.")

    # --- Print Final Alignment for User Confidence ---
    print("\n>>> Final Sample Integrity Check (First 5 Samples):")
    final_geno = get_samples(args.genotype_matrix, sep)
    final_expr = get_samples(args.gene_expression_matrix, sep)
    
    print(f"  Genotype  ({len(final_geno)}): {final_geno[:5]} ...")
    print(f"  Expression({len(final_expr)}): {final_expr[:5]} ...")
    
    if args.covariates:
        final_cov = get_samples(args.covariates, sep)
        print(f"  Covariates({len(final_cov)}): {final_cov[:5]} ...")
    else:
        print("  Covariates: None")
    print("-" * 50 + "\n")


def main():
    args = get_args()

    # Auto-align inputs to avoid null results
    align_inputs(args)

    # Fix QQ Plot Location
    # If using default filename, place it in the output directory
    if args.qq_plot == 'MatrixEqtlQQPlot.pdf':
        out_dir = os.path.dirname(args.output_file)
        if not out_dir: out_dir = '.'
        args.qq_plot = os.path.join(out_dir, "qq_plot.pdf")
    
    print(f">>> Full Analysis Start.")
    print(f"    Results will be saved to: {args.output_file}")
    print(f"    QQ Plot will be saved to: {args.qq_plot}")

    r = ro.r

    # Determine directory where this script resides
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Check Mode
    has_loc = (args.genotype_positions is not None) and (args.gene_positions is not None)

    if has_loc:
        print(">>> Running Standard Mode (Cis/Trans with Locations)...")
        r.source(os.path.join(script_dir, "mxeqtl.R"))
        r.mxeqtl(args.genotype_matrix, args.genotype_positions, args.gene_expression_matrix, args.gene_positions,
                 covariates=args.covariates, cis_output_file=args.output_file, cis_pval=args.p_value, trans_output_file=args.trans_output_file,
                 trans_pval=args.trans_p_value, cis_dist=args.cis_distance, MAF=args.maf, qq=args.qq_plot, model=args.model,
                 header=args.no_header, rownames=args.no_rownames, missing=args.missing, sep=args.sep, chunk_size=args.chunk_size,
                 no_fdr=args.no_fdr)
    
    else:
        print(">>> Running Location-Free Mode (Classic All-vs-All)...")
        print("    Note: Cis-distance ignored. Treating all pairs as Trans.")
        r.source(os.path.join(script_dir, "mxeqtl_noloc.R"))
        
        # Determine strict output file/pval (Prefer Trans since logic is Trans)
        target_out = args.trans_output_file if args.trans_output_file else args.output_file
        target_pval = args.trans_p_value if args.trans_p_value > 0 else args.p_value
        
        r.run_engine(
            snp_file = args.genotype_matrix,
            expr_file = args.gene_expression_matrix,
            cov_file = args.covariates,
            output_file = target_out,
            p_value = target_pval,
            chunk_size = args.chunk_size,
            no_fdr = args.no_fdr
        )

if __name__ == '__main__':
    main()
