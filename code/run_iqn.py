from __future__ import print_function
import argparse
import os
import sys
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame
import numpy
from check_file import check_file


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', required=True, help='')
    parser.add_argument('-o', '--output-file', help='')
    parser.add_argument('-s', '--stdout', action='store_true', help='')
    args = parser.parse_args()

    check_file(args.input_file)

    if not args.stdout and not args.output_file:
        args.output_file = args.input_file + '.qnorm'

    return(args.input_file, args.output_file)


def iqn(input_filename, output_filename=None):
    r = ro.r
    utils = importr('utils', robject_translations={'with': '_with'})
    write_table = utils.write_table

    
    # Determine directory where this script resides
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Source R scripts relative to this directory
    r.source(os.path.join(script_dir, "general_iqn_py.R"))

    if output_filename:
        out_dir = os.path.dirname(output_filename)
        # print(f"DEBUG: Output Dir: {out_dir}") 
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

    print("DEBUG: Running inverse_quantile_norm in R...")   
    normed = r['inverse_quantile_norm'](input_filename)
    
    print(f"DEBUG: Writing table to {output_filename}...")
    write_table(normed, output_filename, col_names=True, row_names=False, quote=False, sep="\t")
    
    if os.path.exists(output_filename):
        print(f"DEBUG: SUCCESS - File verified at {output_filename}")
        print(f"DEBUG: File size: {os.path.getsize(output_filename)} bytes")
        
        # --- Automatic Verification & Plotting ---
        try:
            import matplotlib
            matplotlib.use('Agg') # Headless mode
            import matplotlib.pyplot as plt
            import pandas as pd
            import numpy as np

            print("DEBUG: Generating normality check plot...")
            # Read first 1000 rows for quick check
            df = pd.read_csv(output_filename, sep='\t', index_col=0, nrows=1000)
            
            # Plot
            plot_file = output_filename + ".check.png"
            plt.figure(figsize=(10, 6))
            
            # Plot distribution of 3 random genes
            if len(df) > 3:
                sample_genes = df.sample(3).index
            else:
                sample_genes = df.index
            
            for gene in sample_genes:
                data = df.loc[gene]
                plt.hist(data, bins=20, alpha=0.5, density=True, label=f'Gene: {gene}')
            
            # Standard Normal overlay
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 100)
            p = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x**2)
            plt.plot(x, p, 'k--', linewidth=2, label='Standard Normal (0,1)')
            
            plt.title('IQN Result Verification (Random Genes)')
            plt.xlabel('Normalized Expression Value')
            plt.ylabel('Density')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(plot_file)
            print(f"DEBUG: Normality plot saved to {plot_file}")
            
        except ImportError:
            print("WARNING: matplotlib or pandas not found. Skipping plot generation.")
        except Exception as e:
            print(f"WARNING: Failed to generate plot: {e}")
        # -----------------------------------------

    else:
        print(f"DEBUG: FAILURE - File NOT found at {output_filename}")


def main():
    input_filename, output_filename = get_args()
    iqn(input_filename, output_filename)

if __name__ == '__main__':
    main()
