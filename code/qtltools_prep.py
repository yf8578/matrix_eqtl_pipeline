import pandas as pd
import argparse
import sys
import os
import subprocess

def check_dependencies():
    """Check if bgzip and tabix are available."""
    for cmd in ['bgzip', 'tabix']:
        if subprocess.call(['which', cmd], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"Error: {cmd} is not found in PATH. Please install htslib/tabix.")
            sys.exit(1)

def get_args():
    parser = argparse.ArgumentParser(description="Prepare Expression Data for QTLtools (BED format + Tabix Index)")
    parser.add_argument('--expression', required=True, help='Expression Matrix file (ID in first col)')
    parser.add_argument('--positions', help='Gene Position file (geneid, chr, left, right). Optional if using --dummy-pos')
    parser.add_argument('--out', required=True, help='Output BED file prefix (e.g. output/qtltools_input/expr)')
    parser.add_argument('--dummy-pos', action='store_true', help='Generate dummy positions (for simple traits/metabolites)')
    parser.add_argument('--dummy-chr', default='1', help='Chromosome for dummy positions')
    return parser.parse_args()

def main():
    args = get_args()
    check_dependencies()
    
    # Create output dir if needed
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    final_bed = args.out
    if not final_bed.endswith('.bed'):
        final_bed += '.bed'

    print(f"Reading expression: {args.expression}")
    expr = pd.read_csv(args.expression, sep='\t')
    # Use the first column as ID (rename to geneid for consistency)
    id_col = expr.columns[0]
    expr.rename(columns={id_col: 'geneid'}, inplace=True)

    # DataFrame to construct BED
    # Format: #Chr, start, end, pid, gid, strand, samples...
    
    if args.dummy_pos:
        print(f"Generating dummy positions on chr{args.dummy_chr}...")
        n_feat = len(expr)
        # 1000, 2000, 3000...
        starts = [1000 + i*1000 for i in range(n_feat)]
        
        bed = pd.DataFrame()
        bed['#Chr'] = [args.dummy_chr] * n_feat
        bed['start'] = starts
        bed['end'] = [s + 100 for s in starts]
        bed['pid'] = expr['geneid']
        bed['gid'] = expr['geneid']
        bed['strand'] = '.'
        
        # Add sample data
        bed = pd.concat([bed, expr.iloc[:, 1:]], axis=1)
        
    else:
        if not args.positions:
            print("Error: --positions file is required unless --dummy-pos is specified.")
            sys.exit(1)
            
        print(f"Reading positions: {args.positions}")
        pos = pd.read_csv(args.positions, sep='\t')
        # Normalize headers
        pos.columns = [c.lower() for c in pos.columns]
        if 'features' in pos.columns: pos.rename(columns={'features': 'geneid'}, inplace=True)
        
        # Merge
        print("Merging expression and positions...")
        merged = pd.merge(pos, expr, on='geneid', how='inner')
        print(f"  Matched {len(merged)} genes.")
        
        bed = pd.DataFrame()
        bed['#Chr'] = merged['chr'].astype(str)
        bed['start'] = merged['left'].astype(int)
        bed['end'] = merged['right'].astype(int)
        bed['pid'] = merged['geneid']
        bed['gid'] = merged['geneid']
        bed['strand'] = '.'
        
        # Add sample columns
        sample_cols = [c for c in expr.columns if c != 'geneid']
        bed = pd.concat([bed, merged[sample_cols]], axis=1)

    # Sort
    print("Sorting BED file...")
    # Sort key: Chromosome (str) then Start (int)
    # Note: Traditional sort -k1,1 -k2,2n is handled by pandas sort_values
    bed.sort_values(by=['#Chr', 'start'], ascending=[True, True], inplace=True)

    # Save Uncompressed
    print(f"Writing to {final_bed}...")
    bed.to_csv(final_bed, sep='\t', index=False, header=True, float_format='%.5g')

    # bgzip
    print("Compressing with bgzip...")
    subprocess.check_call(f"bgzip -f {final_bed}", shell=True)
    gz_file = final_bed + ".gz"
    
    # tabix
    print("Indexing with tabix...")
    subprocess.check_call(f"tabix -p bed {gz_file}", shell=True)

    print(f"Success! Output ready: {gz_file}")

if __name__ == "__main__":
    main()
