import pandas as pd
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Check normality of expression matrix")
    parser.add_argument('-i', '--input', required=True, help='Path to expression matrix (qnorm)')
    parser.add_argument('-o', '--plot-output', default='normality_check.png', help='Output path for histogram plot')
    return parser.parse_args()

def main():
    args = get_args()
    
    print(f"Reading {args.input}...")
    try:
        # Read first 1000 rows to save time, or full file if not too big. 
        # Expression matrices are usually wide (samples) x tall (genes).
        # IQN normalizes *each gene* (row) across samples.
        df = pd.read_csv(args.input, sep='\t', index_col=0, nrows=1000)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    print(f"Dataset shape: {df.shape}")
    
    # Check stats for the first few genes
    means = df.mean(axis=1)
    stds = df.std(axis=1)
    
    print("\n--- Summary Statistics (First 5 Genes) ---")
    print(f"{'GeneID':<20} {'Mean':<10} {'Std Dev':<10}")
    for i in range(min(5, len(df))):
        print(f"{df.index[i]:<20} {means.iloc[i]:.4f}     {stds.iloc[i]:.4f}")
        
    print("\n--- Overall Distribution Check ---")
    print(f"Average Mean across all loaded genes: {means.mean():.4f} (Expected ~0)")
    print(f"Average Std Dev across all loaded genes: {stds.mean():.4f} (Expected ~1)")

    # Plot histogram for a random gene
    if args.plot_output:
        plt.figure(figsize=(10, 6))
        
        # Pick 3 random genes
        if len(df) > 3:
            sample_genes = df.sample(3).index
        else:
            sample_genes = df.index
            
        for gene in sample_genes:
            data = df.loc[gene]
            plt.hist(data, bins=20, alpha=0.5, label=str(gene), density=True)
            
        # Plot standard normal for comparison
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x**2)
        plt.plot(x, p, 'k', linewidth=2, label='Standard Normal')
        
        plt.title('Distribution of Normalized Expression (Random 3 Genes)')
        plt.xlabel('Value')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(args.plot_output)
        print(f"\nSaved distribution plot to: {args.plot_output}")

if __name__ == '__main__':
    main()
