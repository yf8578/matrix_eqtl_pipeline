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

    if output_filename is None:
        output_filename = ""

    if output_filename:
        out_dir = os.path.dirname(output_filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

    normed = r['inverse_quantile_norm'](input_filename)
    
    write_table(normed, output_filename, col_names=True, row_names=False, quote=False, sep="\t")


def main():
    input_filename, output_filename = get_args()
    iqn(input_filename, output_filename)

if __name__ == '__main__':
    main()
