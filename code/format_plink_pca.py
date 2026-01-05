import pandas as pd
import argparse
import os
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Convert PLINK .eigenvec to MatrixEQTL covariate format.")
    parser.add_argument('-i', '--input', required=True, help='Path to PLINK .eigenvec file')
    parser.add_argument('-o', '--output', required=True, help='Path to output file (MatrixEQTL format)')
    parser.add_argument('-n', '--num-pcs', type=int, default=3, help='Number of PCs to extract (default: 3)')
    return parser.parse_args()

def main():
    args = get_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found.")
        sys.exit(1)

    try:
        # Read PLINK .eigenvec (space-separated, no header)
        # Expected columns: FID, IID, PC1, PC2, ...
        # If headers exist (PLINK 2.0+), pandas usually handles it if we account for it, 
        # but standard PLINK 1.9 .eigenvec is usually clean. 
        # Using header=None safely.
        df = pd.read_csv(args.input, sep='\s+', header=None)
        
        # Verify shape
        if df.shape[1] < 2 + args.num_pcs:
            print(f"Error: Input file has {df.shape[1]} columns, expected at least {2 + args.num_pcs} (2 ID cols + {args.num_pcs} PCs).")
            sys.exit(1)

        # Extract SampleID (col 1, index 1) and PCs (cols 2 to 2+N)
        # Col 0 is FID, Col 1 is IID.
        # Data starts from col 2.
        df_out = df.iloc[:, 2 : 2 + args.num_pcs]
        
        # Create column names PC1, PC2...
        df_out.columns = [f"PC{i+1}" for i in range(args.num_pcs)]
        
        # Set SampleID as index
        df_out.index = df.iloc[:, 1]
        
        # Transpose: MatrixEQTL expects Row=Covariate, Col=Sample
        df_transposed = df_out.transpose()
        
        # Ensure output dir exists
        out_dir = os.path.dirname(args.output)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        # Save
        df_transposed.to_csv(args.output, sep='\t', index=True, index_label='ID')
        print(f"Successfully converted PCA. Size: {df_transposed.shape}. Saved to: {args.output}")
        
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
