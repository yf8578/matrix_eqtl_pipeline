import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

def combine_and_plot(pca_file, peer_file, known_file, out_dir, annot=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    dfs = []
    labels = []
    
    # helper to read
    def read_cov(fpath, label, transpose_if_needed=False, samples_for_check=None):
        if not fpath or not os.path.exists(fpath):
            if fpath: print(f"Warning: File not found for {label}: {fpath}")
            return None
            
        print(f"Loading {label}: {fpath}")
        
        # Try sniffing format (CSV with spaces vs TSV)
        # Using python engine for robust regex separator if needed, or simple try/except
        try:
            # First try tab
            df = pd.read_csv(fpath, sep='\t', index_col=0)
            if df.shape[1] <= 1: # Suspicion: might be CSV?
                # Reset and try comma
                df_csv = pd.read_csv(fpath, sep=',', index_col=0)
                if df_csv.shape[1] > df.shape[1]:
                    df = df_csv
                    print(f"  (Detected CSV format)")
        except:
            # Fallback to defaults
            df = pd.read_csv(fpath, sep=None, engine='python', index_col=0)

        # Cleanup whitespace in Index and Columns if they are strings
        if df.index.dtype == 'object':
            df.index = df.index.str.strip()
        if df.columns.dtype == 'object':
            df.columns = df.columns.str.strip()
            
        # Drop empty index rows (like the ones in user screenshot)
        df = df[df.index.notna() & (df.index != '')]
        
        # Transpose Check
        # MatrixEQTL needs: Rows=Covariates, Cols=Samples
        # If transpose_if_needed is True, we check if the Index looks like Samples (present in other files)
        if transpose_if_needed and samples_for_check is not None:
            # Check overlap of Index with Reference Samples
            index_intersect = len(set(df.index).intersection(samples_for_check))
            col_intersect = len(set(df.columns).intersection(samples_for_check))
            
            if index_intersect > col_intersect:
                print(f"  Detected Sample x Feature orientation. Transposing '{label}'.")
                df = df.T
        
        # Deduplicate Columns (Samples)
        if not df.columns.is_unique:
            print(f"  Warning: Duplicate columns (samples) found in '{label}'. Keeping first occurrence.")
            df = df.loc[:, ~df.columns.duplicated()]
            
        # Deduplicate Rows (Covariates/Features) - Just in case
        if not df.index.is_unique:
             print(f"  Warning: Duplicate rows found in '{label}'. Keeping first occurrence.")
             df = df.loc[~df.index.duplicated(), :]
             
        # Drop columns with empty names
        cols_to_drop = [c for c in df.columns if str(c).strip() == '' or str(c).lower() == 'nan']
        if cols_to_drop:
            print(f"  Dropping {len(cols_to_drop)} columns with empty/invalid names.")
            df = df.drop(columns=cols_to_drop)
            
        # FORCE NUMERIC: Coerce anything that looks like a number
        # (Handles strings like " 1", "2 " etc which might have slipped in)
        df = df.apply(pd.to_numeric, errors='coerce')
        
        # Drop Constant Rows (Variance = 0) - e.g. Disease=1,1,1...
        # Correlation is undefined for constant variables.
        std = df.std(axis=1) # std across samples
        constant_features = std[std == 0].index
        if len(constant_features) > 0:
            print(f"  Warning: Dropping {len(constant_features)} constant features (Variance=0): {list(constant_features)}")
            print("           (Constant covariates cannot be used for correlation/analysis)")
            df = df.drop(index=constant_features)
            
        print(f"  Shape: {df.shape}")
        return df

    # Load Base (PCA/PEER) to get valid samples list
    df_pca = read_cov(pca_file, "Genotype PCs")
    valid_samples = set(df_pca.columns) if df_pca is not None else None
    
    df_peer = read_cov(peer_file, "PEER Factors")
    if valid_samples is None and df_peer is not None:
         valid_samples = set(df_peer.columns)
    
    df_known = read_cov(known_file, "Known Covariates", transpose_if_needed=True, samples_for_check=valid_samples)
    
    # Append
    if df_pca is not None: dfs.append(df_pca)
    if df_peer is not None: dfs.append(df_peer)
    if df_known is not None: dfs.append(df_known)
    
    if not dfs:
        print("Error: No valid covariate files provided.")
        sys.exit(1)
        
    # Merge (Inner Join on Samples aka Columns)
    # MatrixEQTL format: Rows = Covariates, Cols = Samples
    # So we align columns.
    
    final_df = pd.concat(dfs, axis=0, join='inner')
    
    print(f"Combined Covariates Shape: {final_df.shape}")
    print(f"Samples: {final_df.shape[1]}, Covariates: {final_df.shape[0]}")
    
    # Save Combined
    # Ensure index has a name so to_csv writes a header for it (Crucial for QTLtools)
    final_df.index.name = 'id'
    
    out_file = out_dir / "final_covariates.txt"
    final_df.to_csv(out_file, sep='\t')
    print(f"Saved combined covariates to: {out_file}")
    
    # Correlation Analysis
    # Corr needs Samples as rows
    corr_matrix = final_df.T.corr()
    
    # Save Corr Matrix
    corr_file = out_dir / "covariates_correlation.txt"
    corr_matrix.to_csv(corr_file, sep='\t')
    
    # Plotting
    try:
        # Dynamic figsize based on number of variables
        n_vars = len(corr_matrix.columns)
        # Estimate size: at least 10x10, but grow with n_vars
        plot_size = max(12, n_vars * 0.3)
        
        plt.figure(figsize=(plot_size, plot_size * 0.8))
        sns.heatmap(corr_matrix, 
                    annot=annot, 
                    cmap='RdBu_r', 
                    center=0, 
                    square=True,
                    fmt=".2f",
                    linewidths=.5 if n_vars < 50 else 0,
                    cbar_kws={"shrink": .8},
                    xticklabels=True, # Force showing all labels
                    yticklabels=True) # Force showing all labels
                    
        plt.title("Covariate Correlation Matrix")
        plt.tight_layout()
        
        plot_file = out_dir / "covariates_correlation.pdf"
        plt.savefig(plot_file)
        print(f"Saved correlation plot to: {plot_file}")
        
    except Exception as e:
        print(f"Warning: Plotting failed (maybe missing libraries?): {e}")

def main():
    parser = argparse.ArgumentParser(description="Combine Covariates & Plot Correlation")
    parser.add_argument("--pca", help="Genotype PCs file")
    parser.add_argument("--peer", help="PEER Factors file")
    parser.add_argument("--known", help="Known Covariates file")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--no-annot", action='store_true', help="Do not show correlation values on the heatmap")
    
    args = parser.parse_args()
    
    combine_and_plot(args.pca, args.peer, args.known, args.out_dir, annot=(not args.no_annot))

if __name__ == "__main__":
    main()
