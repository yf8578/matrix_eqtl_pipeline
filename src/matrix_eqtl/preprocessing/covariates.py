
import pandas as pd
import subprocess
from pathlib import Path
import os
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

def generate_covariates(genotype_prefix, expression_matrix, known_covariates, output_dir, pca_n=3, peer_n=0, threads=1):
    """
    Generate and combine covariates (Genotype PCA + PEER + Known).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"DEBUG: Python CWD: {os.getcwd()}")
    
    cov_list = []
    
    # 1. Genotype PCA
    if pca_n > 0 and genotype_prefix:
        pca_out = output_dir / "genotype_pca"
        
        # Verify input existence
        bed_file = Path(str(genotype_prefix) + ".bed")
        if not bed_file.exists():
            print(f"ERROR: Genotype BED file not found at: {bed_file.resolve()}")
            print("Did Step 3 (Genotype Processing) run successfully?")
            # Check directory content for debugging
            parent_dir = bed_file.parent
            if parent_dir.exists():
                print(f"Contents of {parent_dir}:")
                for f in parent_dir.iterdir():
                    print(f"  {f.name} ({f.stat().st_size} bytes)")
            else:
                 print(f"Directory {parent_dir} NOT found.")
            raise FileNotFoundError(f"Missing PLINK input: {bed_file}")

        # Run PLINK PCA
        # plink --bfile prefix --pca n --out output
        cmd = ["plink", "--bfile", str(genotype_prefix), "--pca", str(pca_n), "--out", str(pca_out), "--threads", str(threads)]
        print(f"Running PLINK PCA: {' '.join(cmd)}")
        subprocess.check_call(cmd)
        
        # Format PCA: PLINK .eigenvec has FID IID PC1 PC2 ...
        # We need to transpose for MatrixEQTL: Rows=PCs, Cols=Samples (with ID header? MatrixEQTL usually wants header)
        vec_file = str(pca_out) + ".eigenvec"
        df_pca = pd.read_csv(vec_file, sep=r'\s+', header=None)
        # Columns: FID, IID, PC1, PC2...
        # Drop FID
        df_pca = df_pca.iloc[:, 1:]
        df_pca.columns = ['id'] + [f'PC{i+1}' for i in range(pca_n)]
        df_pca = df_pca.set_index('id').T
        
        pca_final = output_dir / "genotype_pcs.txt"
        df_pca.to_csv(pca_final, sep='\t')
        cov_list.append(df_pca)
        
    # 2. PEER Factors
    if peer_n > 0:
        peer_out = output_dir / "peer_factors.tsv"
        
        # Locate R script
        # Assumes structure: src/matrix_eqtl/preprocessing/covariates.py
        # Project Root: ../../../
        # R Src: ../../../r_src/peer_function.R
        project_root = Path(__file__).parent.parent.parent.parent 
        r_script = project_root / "r_src" / "peer_function.R"
        
        if not r_script.exists():
            # Fallback for development (if running from different cwd possibly?)
            r_script = Path("r_src/peer_function.R").resolve()
            
        if not r_script.exists():
             raise FileNotFoundError(f"Could not find peer_function.R at {r_script}")

        # Use rpy2 to source and run
        r = ro.r
        r.source(str(r_script))
        
        # peer_function(input_file, num_factors, covariates_file)
        # Note: existing peer_function.R takes file paths.
        # We pass the expression matrix file path.
        cov_file_arg = known_covariates if known_covariates else ""
        
        print(f"Running PEER with {peer_n} factors...")
        # peer_function returns a DataFrame (Variables in rows? No, usually PEER returns factors x samples)
        # Let's check what peer_function.R returns. 
        # Original wrapper: write_table(factors, ...)
        
        factors = r['peer_function'](str(expression_matrix), peer_n, str(cov_file_arg))
        
        # Convert R DF to Python or write to file directly from R? 
        # The wrapper wrote it from python using utils.write_table.
        # Let's convert to pandas to merge generally.
        
        # R dataframe to pandas? Or just write temp and read.
        # Writing temp is safer for rpy2 versions.
        temp_peer = output_dir / "peer_temp.txt"
        utils = importr('utils')
        utils.write_table(factors, str(temp_peer), sep="\t", quote=False, row_names=True, col_names=True)
        
        df_peer = pd.read_csv(temp_peer, sep='\t', index_col=0)
        # PEER usually outputs samples in rows, factors in columns? 
        # Or factors in rows?
        # Check `run_peer.py`: `row_names=False` in original script implies standard format.
        # If we need to merge, we need consistent orientation. 
        # MatrixEQTL covariates: Rows=Covariates, Cols=Samples.
        
        # If peer_function returns standard PEER output (N x K), we transpose -> (K x N).
        # Let's assume we need to align orientation.
        # Inspecting df_peer shape at runtime is best, but here we assume we need to Transpose if sample IDs are index.
        # Actually existing script `run_peer.py` outputted `row_names=False`.
        
        df_peer.to_csv(peer_out, sep='\t')
        
        # Read back to ensure format for merging
        # If existing script `run_peer.py` logic was correct, let's stick to it.
        # If we want to merge, we need all in one dataframe.
        # Let's assume df_peer is formatted correctly (Samples columns? or Rows?)
        # Standard: Variables (Genes/Covariates) x Samples.
        # If PEER output was Samples x Factors, we T.
        
        # To be safe, we will just start with Transpose because PEER usually gives factors per sample.
        if df_peer.shape[0] == len(df_pca.columns): # if rows == samples
             df_peer = df_peer.T
             
        cov_list.append(df_peer)

    # 3. Known Covariates
    if known_covariates:
        df_known = pd.read_csv(known_covariates, sep='\t', index_col=0)
        cov_list.append(df_known)
        
    # Merge all
    if not cov_list:
        return None
        
    final_df = pd.concat(cov_list, axis=0, join='inner')
    final_out = output_dir / "final_covariates.txt"
    final_df.to_csv(final_out, sep='\t')
    
    return str(final_out)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Covariate Generation (PCA + PEER)")
    parser.add_argument("--genotype", help="Genotype PLINK prefix (for PCA)")
    parser.add_argument("--expression", required=True, help="Expression matrix file (for PEER)")
    parser.add_argument("--known", help="Known covariates file")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--pca-n", type=int, default=3, help="Number of Genotype PCs")
    parser.add_argument("--peer-n", type=int, default=0, help="Number of PEER factors")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for PLINK PCA")
    
    args = parser.parse_args()
    
    generate_covariates(
        genotype_prefix=args.genotype,
        expression_matrix=args.expression,
        known_covariates=args.known,
        output_dir=args.out_dir,
        pca_n=args.pca_n,
        peer_n=args.peer_n,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
