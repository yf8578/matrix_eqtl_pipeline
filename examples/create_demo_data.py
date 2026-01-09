import pandas as pd
import argparse
import os
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Create Demo Dataset by Subsampling")
    parser.add_argument('--genotype', required=True)
    parser.add_argument('--snp-pos', required=True)
    parser.add_argument('--expression', required=True)
    parser.add_argument('--gene-pos', required=True)
    parser.add_argument('--covariates', required=True)
    
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--n-samples', type=int, default=50)
    parser.add_argument('--n-snps', type=int, default=2000)
    parser.add_argument('--n-genes', type=int, default=100)
    parser.add_argument('--seed', type=int, default=42)
    
    return parser.parse_args()

def subset_matrix(fpath, samples, n_rows, seed, name):
    print(f"  Subsetting {name} ({n_rows} rows)...")
    df = pd.read_csv(fpath, sep='\t')
    
    # 1. Filter Samples
    id_col = df.columns[0]
    # Intersection of requested samples and present samples
    present = [s for s in samples if s in df.columns]
    
    if len(present) < len(samples):
        print(f"    Warning: Only {len(present)}/{len(samples)} requested samples found in {name}.")
    
    # Random row selection
    if len(df) > n_rows:
        df_sub = df.sample(n=n_rows, random_state=seed)
    else:
        df_sub = df
        
    # Select cols
    cols = [id_col] + present
    return df_sub[cols]

def subset_pos(fpath, ids, name):
    print(f"  Subsetting {name} positions...")
    df = pd.read_csv(fpath, sep='\t')
    # Filter by IDs (first column usually)
    id_col = df.columns[0]
    # Ensure IDs are strings
    df[id_col] = df[id_col].astype(str)
    
    sub = df[df[id_col].isin(ids)]
    return sub

def main():
    args = get_args()
    if not os.path.exists(args.out_dir): os.makedirs(args.out_dir)
    
    print(">>> Creating Demo Data...")
    
    # 1. Pick Samples from Genotype (Master)
    print("  Reading Master Genotype Header...")
    geno_header = pd.read_csv(args.genotype, sep='\t', nrows=1).columns.tolist()
    all_samples = geno_header[1:]
    
    # Sample subset
    import random
    random.seed(args.seed)
    if len(all_samples) > args.n_samples:
        selected_samples = random.sample(all_samples, args.n_samples)
    else:
        selected_samples = all_samples
    
    print(f"  Selected {len(selected_samples)} samples.")
    
    # 2. Subset Genotype Matrix
    geno_sub = subset_matrix(args.genotype, selected_samples, args.n_snps, args.seed, "Genotype")
    geno_sub.to_csv(os.path.join(args.out_dir, "genotype.demo.txt"), sep='\t', index=False)
    
    # 3. Subset Genotype Pos
    # IDs in geno_sub
    snp_ids = geno_sub.iloc[:, 0].astype(str).tolist()
    snp_pos_sub = subset_pos(args.snp_pos, snp_ids, "SNP Pos")
    snp_pos_sub.to_csv(os.path.join(args.out_dir, "snp_pos.demo.txt"), sep='\t', index=False)
    
    # 4. Subset Expression Matrix
    expr_sub = subset_matrix(args.expression, selected_samples, args.n_genes, args.seed, "Expression")
    expr_sub.to_csv(os.path.join(args.out_dir, "expression.demo.txt"), sep='\t', index=False)
    
    # 5. Subset Gene Pos
    gene_ids = expr_sub.iloc[:, 0].astype(str).tolist()
    gene_pos_sub = subset_pos(args.gene_pos, gene_ids, "Gene Pos")
    gene_pos_sub.to_csv(os.path.join(args.out_dir, "gene_pos.demo.txt"), sep='\t', index=False)
    
    # 6. Subset Covariates
    # Covariates usually few rows, keep all rows, subset columns
    print("  Subsetting Covariates (Samples only)...")
    cov = pd.read_csv(args.covariates, sep='\t')
    cov_id = cov.columns[0]
    cov_present = [s for s in selected_samples if s in cov.columns]
    cov_sub = cov[[cov_id] + cov_present]
    cov_sub.to_csv(os.path.join(args.out_dir, "covariates.demo.txt"), sep='\t', index=False)
    
    print(f"\n>>> Done! Demo data in: {args.out_dir}")

if __name__ == "__main__":
    main()
