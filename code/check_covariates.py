import pandas as pd
import numpy as np
import argparse
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

def get_args():
    parser = argparse.ArgumentParser(description="Check correlations between Covariates, PEER factors, and Expression PCs")
    parser.add_argument('--expression', required=True, help='Expression matrix (qnorm)')
    parser.add_argument('--peer', required=True, help='PEER factors file')
    parser.add_argument('--covariates', required=True, help='Known covariates file (e.g. Sex, Age)')
    parser.add_argument('--output', default='covariate_check', help='Output prefix')
    return parser.parse_args()

def main():
    args = get_args()
    
    print("Loading data...")
    try:
        # Load Expression (Genes x Samples) -> Transpose to (Samples x Genes) for PCA
        expr = pd.read_csv(args.expression, sep='\t', index_col=0).transpose()
        
        # Load PEER (Factors x Samples) -> Transpose to (Samples x Factors)
        # Note: PEER output usually has 'ID' column as first column. check_file format?
        # Our pipeline output has ID as col 1.
        peer = pd.read_csv(args.peer, sep='\t', index_col=0).transpose()
        
        # Load Known Covariates (Covariates x Samples) -> Transpose
        covs = pd.read_csv(args.covariates, sep='\t', index_col=0).transpose()
        
        # Align Samples
        # Find common samples
        common = expr.index.intersection(peer.index).intersection(covs.index)
        print(f"Common samples: {len(common)}")
        
        if len(common) < 10:
            print("Error: Too few overlapping samples. Check sample IDs.")
            sys.exit(1)
            
        expr = expr.loc[common]
        peer = peer.loc[common]
        covs = covs.loc[common]
        
    except Exception as e:
        print(f"Error loading files: {e}")
        sys.exit(1)

    # 1. Calculate Expression PCs (Top 5)
    print("Computing Expression PCs...")
    pca = PCA(n_components=5)
    pcs = pca.fit_transform(expr)
    pc_df = pd.DataFrame(pcs, index=common, columns=[f"Exp_PC{i+1}" for i in range(5)])
    
    # 2. Combine all factors for correlation
    # We want to see: Covariates vs PEERs, Covariates vs Exp_PCs
    # Let's combine them all into one big dataframe
    # Rename PEER cols to avoid confusion if needed (they are likely PEER1, PEER2...)
    
    combined = pd.concat([covs, peer, pc_df], axis=1)
    
    # Force all to numeric (coerce errors)
    combined = combined.apply(pd.to_numeric, errors='coerce')
    combined = combined.dropna(axis=1, how='all') # Drop non-numeric cols
    
    # 3. Correlation Matrix
    corr_matrix = combined.corr(method='pearson')
    
    # We are mainly interested in: Rows=Known Covariates, Cols=PEERs + PCs
    # But full heatmap is also informative.
    
    # 4. Plot Heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, cmap='coolwarm', center=0, annot=False)
    plt.title("Correlation Heatmap: Known Covariates vs Inferred Factors")
    
    out_plot = args.output + ".heatmap.png"
    plt.tight_layout()
    plt.savefig(out_plot)
    print(f"Saved correlation heatmap: {out_plot}")
    
    # 5. Text Report: High Correlations with Known Covariates
    out_txt = args.output + ".report.txt"
    with open(out_txt, 'w') as f:
        f.write("=== Covariate Correlation Report ===\n\n")
        
        # Check each known covariate
        for cov_name in covs.columns:
            if cov_name not in corr_matrix.index: continue
            
            f.write(f"Covariate: {cov_name}\n")
            # Get correlations with Peers and PCs
            corrs = corr_matrix.loc[cov_name]
            
            # Sort by absolute value
            sorted_corrs = corrs.abs().sort_values(ascending=False)
            
            # Print top 5 correlations (excluding itself)
            count = 0
            for name, val in sorted_corrs.items():
                if name == cov_name: continue
                real_val = corrs[name] # get signed value
                if count < 5:
                    f.write(f"  vs {name:<15}: {real_val:.4f}\n")
                    count += 1
            f.write("\n")
            
    print(f"Saved textual report: {out_txt}")
    
    # Print a preview to stdout
    print("\nTop Correlations for Known Covariates:")
    with open(out_txt, 'r') as f:
        print(f.read())

if __name__ == '__main__':
    main()
