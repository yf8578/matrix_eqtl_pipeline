from __future__ import print_function
import argparse
import os
import sys
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame
from check_file import check_file


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', required=True, help='Input expression file')
    parser.add_argument('-n', '--number-factors', required=True, help='Number of PEER factors')
    parser.add_argument('-o', '--output-file', help='Output file path')
    parser.add_argument('-c', '--covariates', help='File with known covariates (optional)')
    args = parser.parse_args()

    check_file(args.input_file)
    if args.covariates:
        check_file(args.covariates)

    args.output_file = args.output_file or (
        args.input_file + '.peer_factors_' + args.number_factors)

    return(args.input_file, args.number_factors, args.output_file, args.covariates)


def main():
    r = ro.r
    utils = importr('utils', robject_translations={'with': '_with'})
    write_table = utils.write_table

    # Determine directory where this script resides
    script_dir = os.path.dirname(os.path.realpath(__file__))

    r.source(os.path.join(script_dir, "peer_function.R"))

    input_filename, num_factors, output_filename, covariates_file = get_args()
    
    if output_filename:
        out_dir = os.path.dirname(output_filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

    # Pass covariates_file (can be None) to R function
    # Convert None to "" to avoid rpy2 conversion error
    r_cov_file = covariates_file if covariates_file else ""
    factors = r['peer_function'](input_filename, num_factors, r_cov_file)

    write_table(factors, output_filename, col_names=True,
                row_names=False, sep="\t", quote=False)

if __name__ == '__main__':
    main()
