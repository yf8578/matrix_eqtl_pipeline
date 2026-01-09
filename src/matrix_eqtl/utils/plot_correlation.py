import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys

def plot_correlation(input_file, output_prefix=None):
    """
    Plots the correlation heatmap for a MatrixEQTL Covariates file.
    Input format: Rows=Covariates, Cols=Samples (Tab-separated)
    """
    input_path = Path(input_file).resolve()
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)
        
    print(f"Reading covariates from: {input_path}")
    # Read (Rows=Covariates, Cols=Samples)
    df = pd.read_csv(input_path, sep='\t', index_col=0)
    
    # 1. Clean Data
    # Coerce to numeric
    df = df.apply(pd.to_numeric, errors='coerce')
    
    # Drop rows that are all NaN (failed conversion)
    df = df.dropna(how='all', axis=0)
    
    # Drop constant rows (Variance=0)
    std = df.std(axis=1)
    constant = std[std == 0].index
    if len(constant) > 0:
        print(f"Dropping {len(constant)} constant covariates: {list(constant)}")
        df = df.drop(index=constant)

    print(f"Final Covariates Shape: {df.shape}")
    print(f"Variables: {list(df.index)}")
    
    # 2. Calculate Correlation
    # Correlation is between Variables (Rows).
    # Pandas corr() calculates correlation between Columns.
    # So we must Transpose: (Samples x Covariates).
    # Result: (Covariates x Covariates) correlation matrix.
    corr_matrix = df.T.corr()
    
    # Save Matrix
    if output_prefix:
        matrix_out_path = Path(str(output_prefix).replace(".pdf", "").replace(".png", "") + "_matrix.txt")
    else:
        matrix_out_path = input_path.parent / (input_path.stem + "_correlation_matrix.txt")
        
    corr_matrix.to_csv(matrix_out_path, sep='\t')
    print(f"Saved Correlation Matrix to: {matrix_out_path}")
    
    # 3. Plot
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, 
                annot=True,       # Show numbers
                fmt=".2f",        # 2 decimal places
                cmap='RdBu_r',    # Red=Pos, Blue=Neg
                center=0, 
                square=True,
                linewidths=.5,
                cbar_kws={"shrink": .8})
    
    plt.title(f"Covariate Correlation\n({input_path.name})")
    plt.tight_layout()
    
    # Output file
    if output_prefix:
        out_path = Path(output_prefix)
    else:
        out_path = input_path.parent / (input_path.stem + "_correlation.pdf")
        
    plt.savefig(out_path)
    print(f"Saved Correlation Heatmap to: {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot Covariate Correlation Heatmap")
    parser.add_argument("input_file", help="Path to covariates file (Rows=Covariates, Cols=Samples)")
    parser.add_argument("--out", help="Output PDF file path (optional)")
    
    args = parser.parse_args()
    
    plot_correlation(args.input_file, args.out)

if __name__ == "__main__":
    main()
