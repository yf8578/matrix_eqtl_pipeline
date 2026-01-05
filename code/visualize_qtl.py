import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import os
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="Visualize eQTL Results")
    parser.add_argument('--results', required=True, help='MatrixEQTL Result File (cis or trans)')
    parser.add_argument('--snp-pos', required=True, help='SNP Position File (snp, chr, pos)')
    parser.add_argument('--genotype', help='Genotype Matrix (for Boxplots)')
    parser.add_argument('--expression', help='Expression Matrix (for Boxplots)')
    parser.add_argument('--out-prefix', required=True, help='Output prefix (e.g., results/plots/cis)')
    parser.add_argument('--top-n', type=int, default=5, help='Number of top hits to boxplot')
    return parser.parse_args()

def plot_manhattan(results, snp_pos, out_file):
    print("  Drawing Manhattan Plot...")
    
    # Merge Results with Positions
    # Results: SNP, gene, beta, t-stat, p-value, FDR
    # SNP Pos: snp, chr, pos
    
    merged = pd.merge(results, snp_pos, left_on='SNP', right_on='snp', how='inner')
    
    if merged.empty:
        print("    Warning: No SNPs matched positions. Skipping Manhattan.")
        return

    # Prepare data
    merged['p-value'] = pd.to_numeric(merged['p-value'], errors='coerce')
    merged['logp'] = -np.log10(merged['p-value'])
    merged['chr'] = merged['chr'].astype(str)
    
    # Handle Chromosomes (sort numerically X, Y, MT)
    def clean_chr(c):
        c = c.replace('chr', '')
        if c == 'X': return 23
        if c == 'Y': return 24
        if c == 'MT': return 25
        try: return int(c)
        except: return 99
        
    merged['chr_num'] = merged['chr'].apply(clean_chr)
    merged.sort_values(by=['chr_num', 'pos'], inplace=True)
    
    # Plotting
    plt.figure(figsize=(12, 6))
    
    # Color mapping for alt chromosomes
    colors = ['#5D8AA8', '#E32636']
    x_labels = []
    x_labels_pos = []
    
    offset = 0
    for i, c in enumerate(sorted(merged['chr_num'].unique())):
        subset = merged[merged['chr_num'] == c]
        if subset.empty: continue
        
        subset['rel_pos'] = subset['pos'] + offset
        
        plt.scatter(subset['rel_pos'], subset['logp'], c=colors[i % 2], s=10, label=f"Chr {c}")
        
        # Center label
        mid_point = offset + (subset['pos'].max() + subset['pos'].min()) / 2
        x_labels.append(str(c))
        x_labels_pos.append(mid_point)
        
        # Update offset (add buffer)
        offset += subset['pos'].max()
    
    plt.xticks(x_labels_pos, x_labels, fontsize=8)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(P-value)')
    plt.title('Manhattan Plot of eQTLs')
    
    # Significance Line (optional, e.g. Bonferroni)
    # plt.axhline(y=-np.log10(5e-8), color='gray', linestyle='--')
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"    Saved to: {out_file}")

def plot_boxplots(results, geno_file, expr_file, top_n, out_dir):
    print(f"  Drawing Boxplots for Top {top_n} Hits...")
    
    if not geno_file or not expr_file:
         print("    Genotype/Expression files not provided. Skipping Boxplots.")
         return

    # Read matrices (Header only first to check)
    # Just read full for simplicity (assuming manageable size)
    geno = pd.read_csv(geno_file, sep='\t', index_col=0)
    expr = pd.read_csv(expr_file, sep='\t', index_col=0)
    
    # Sort results by p-value
    top_hits = results.sort_values('p-value').head(top_n)
    
    for idx, row in top_hits.iterrows():
        snp = row['SNP']
        gene = row['gene']
        pval = row['p-value']
        
        if snp not in geno.index or gene not in expr.index:
            continue
            
        g_vals = geno.loc[snp]
        e_vals = expr.loc[gene]
        
        # Align samples
        common = g_vals.index.intersection(e_vals.index)
        
        data = pd.DataFrame({
            'Genotype': g_vals[common],
            'Expression': e_vals[common]
        })
        
        # Clean genotype (round to nearest integer 0, 1, 2)
        data['Genotype'] = data['Genotype'].round().astype(int)
        
        plt.figure(figsize=(6, 6))
        sns.boxplot(x='Genotype', y='Expression', data=data, palette='Set2')
        sns.stripplot(x='Genotype', y='Expression', data=data, color='black', alpha=0.5, jitter=0.2)
        
        plt.title(f"{gene} ~ {snp}\nP={pval:.2e}")
        plt.xlabel("Genotype (Dosage)")
        plt.ylabel("Normalized Expression")
        
        fname = f"{out_dir}/boxplot_{gene}_{snp}.pdf"
        plt.savefig(fname)
        plt.close()
        print(f"    Saved: {fname}")

def main():
    args = get_args()
    
    # Create output dir
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    print(">>> Visualizing eQTL Results...")
    
    # Read Results
    res = pd.read_csv(args.results, sep='\t')
    print(f"  Loaded {len(res)} associations.")
    
    # Read SNP Pos (Header: snp, chr, pos)
    # Check headers (MatrixEQTL input format might vary, flexible read)
    spos = pd.read_csv(args.snp_pos, sep='\t')
    # Standardize columns: first is snp, second chr, third pos
    spos.columns = ['snp', 'chr', 'pos'] + list(spos.columns[3:])
    
    # 1. Manhattan
    plot_manhattan(res, spos, args.out_prefix + "_manhattan.png")
    
    # 2. Boxplots
    plot_boxplots(res, args.genotype, args.expression, args.top_n, out_dir)
    
    print(">>> Visualization Done.")

if __name__ == "__main__":
    main()
