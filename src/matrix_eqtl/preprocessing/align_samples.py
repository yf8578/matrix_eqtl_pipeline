import pandas as pd
import argparse
import sys
import os

def get_args():
    parser = argparse.ArgumentParser(description="Align Sample Columns of Expression and Covariates to match Genotype")
    parser.add_argument('--genotype', required=True, help='Genotype Matrix (Master Sample Order)')
    parser.add_argument('--expression', required=True, help='Expression Matrix to align')
    parser.add_argument('--covariates', required=True, help='Covariates Matrix to align')
    parser.add_argument('--out-dir', required=True, help='Output directory for aligned files')
    return parser.parse_args()

def align_file(target_file, master_samples, out_path):
    print(f"Aligning {os.path.basename(target_file)}...")
    
    # Read file
    # Use pandas for handling TSV
    df = pd.read_csv(target_file, sep='\t')
    
    # First column is ID (geneid / covariate_id)
    id_col = df.columns[0]
    
    # Check if all master samples exist in target
    target_samples = df.columns[1:].tolist()
    missing = [s for s in master_samples if s not in target_samples]
    
    if list(master_samples) == target_samples:
        print("  Already aligned. Copying...")
        df.to_csv(out_path, sep='\t', index=False, float_format='%.5g')
        return

    if missing:
        print(f"  Error: {len(missing)} samples from Genotype are missing in {os.path.basename(target_file)}")
        print(f"  Example missing: {missing[:5]}")
        sys.exit(1)
        
    # Reorder
    # New columns = [ID_Col] + [Master_Sample_1, Master_Sample_2, ...]
    new_order = [id_col] + list(master_samples)
    df_aligned = df[new_order]
    
    print(f"  Reordered columns. Writing to {out_path}...")
    df_aligned.to_csv(out_path, sep='\t', index=False, float_format='%.5g')

def main():
    args = get_args()
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
        
    print(">>> Reading Genotype Master Header...")
    # Just read the header to get sample order
    with open(args.genotype, 'r') as f:
        header = f.readline().strip().split('\t')
        # Assuming first col is ID
        master_samples = header[1:]
    
    print(f"  Found {len(master_samples)} samples in Genotype Matrix.")
    print(f"  Order: {master_samples[:5]} ...")
    
    # Define output paths
    out_expr = os.path.join(args.out_dir, "expression.aligned.txt")
    out_cov = os.path.join(args.out_dir, "covariates.aligned.txt")
    
    # Align Expression
    align_file(args.expression, master_samples, out_expr)
    
    # Align Covariates
    align_file(args.covariates, master_samples, out_cov)
    
    print("-" * 30)
    print("✅ Alignment Complete!")
    print(f"New Expression: {out_expr}")
    print(f"New Covariates: {out_cov}")
    print("Please use these files for MatrixEQTL.")

if __name__ == "__main__":
    main()
