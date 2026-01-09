import argparse
import subprocess
from pathlib import Path
import pandas as pd
import os
import sys
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

def run_peer(expression_file, output_dir, peer_n=5, known_covariates=None):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Resolve absolute path for expression file
    expr_path = Path(expression_file).resolve()
    if not expr_path.exists():
        print(f"ERROR: Expression file not found at: {expr_path}")
        print(f"       (Provided input: {expression_file})")
        print(f"       (Current Working Dir: {os.getcwd()})")
        
        # Check if it was a relative path issue
        if not Path(expression_file).is_absolute():
             print("       Hint: Try providing the FULL ABSOLUTE PATH.")
             
        # List dir content if parent exists
        if expr_path.parent.exists():
             print(f"       Contents of {expr_path.parent}:")
             for f in list(expr_path.parent.iterdir())[:5]: # Show first 5
                 print(f"         {f.name}")
        sys.exit(1)
        
    print(f"Resolved Expression File: {expr_path}")
    
    peer_out = output_dir / "peer_factors.tsv"
    
    # Locate R script
    # Path: src/matrix_eqtl/preprocessing/expression_peer.py -> ../../../r_src/peer_function.R
    project_root = Path(__file__).parent.parent.parent.parent 
    r_script = project_root / "r_src" / "peer_function.R"
    
    if not r_script.exists():
        # Fallback
        r_script = Path("r_src/peer_function.R").resolve()
        
    if not r_script.exists():
         print(f"ERROR: Could not find peer_function.R at {r_script}")
         sys.exit(1)

    # Load R script
    print(f"Loading R script: {r_script}")
    r = ro.r
    try:
        r.source(str(r_script))
    except Exception as e:
        print(f"ERROR: Failed to source R script: {e}")
        sys.exit(1)
        
    cov_file_arg = known_covariates if known_covariates else ""
    
    print(f"Running PEER with {peer_n} factors on {expression_file}...")
    
    try:
        # Call R function
        factors = r['peer_function'](str(expression_file), peer_n, str(cov_file_arg))
    except Exception as e:
        print(f"ERROR during PEER execution: {e}")
        sys.exit(1)
        
    # Write output via R or Python
    # Let's write via python to be sure of format
    
    # Convert R DataFrame to Python
    # Using temp file approach for safety with R/Pandas conversion issues
    temp_peer = output_dir / "peer_temp.txt"
    utils = importr('utils')
    utils.write_table(factors, str(temp_peer), sep="\t", quote=False, row_names=True, col_names=True)
    
    df_peer = pd.read_csv(temp_peer, sep='\t', index_col=0)
    
    # R script adds an 'ID' column with PEER1, PEER2...
    # The default read_csv(index_col=0) might read the meaningless row numbers "1","2" as index.
    # If there is an 'ID' column, use it as the index.
    if 'ID' in df_peer.columns:
        df_peer = df_peer.set_index('ID')
        df_peer.index.name = None # Remove the name 'ID' to match other files style
    
    # Transpose if needed
    # PEER usually returns Factors x Samples? Or Samples x Factors?
    # R script says: factors <- t(factors) ... colnames(factors) <- colnames(expr)
    # So R script returns: Rows=Factors, Cols=Samples (ID + Samples)
    # Let's verify shape.
    print(f"PEER Output Shape: {df_peer.shape}")
    
    # MatrixEQTL covariates: Rows=Covariates, Cols=Samples.
    # If the R script already did transpose, we preserve.
    
    df_peer.to_csv(peer_out, sep='\t')
    print(f"PEER Factors saved to: {peer_out}")

def main():
    parser = argparse.ArgumentParser(description="Expression PEER Factors")
    parser.add_argument("--expression", required=True, help="Expression matrix")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--peer-n", type=int, default=5, help="Number of PEER factors")
    parser.add_argument("--known", help="Known covariates file")
    
    args = parser.parse_args()
    
    print(f"DEBUG: Python CWD: {os.getcwd()}")
    
    run_peer(args.expression, args.out_dir, args.peer_n, args.known)

if __name__ == "__main__":
    main()
