
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess

def process_expression(expression_file, samples, output_dir, filter_low_expr=True, normalize_iqn=True, filter_args=None):
    """
    Filter and normalize expression data.
    filter_args: dict with keys 'min_tpm', 'min_samples', 'min_samples_rel', 'max_missing'
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if filter_args is None:
        filter_args = {}
    
    # Defaults
    min_tpm = filter_args.get('min_tpm', 0.1)
    min_samples_abs = filter_args.get('min_samples', 0)
    min_samples_rel = filter_args.get('min_samples_rel', 0.2)
    max_missing = filter_args.get('max_missing', 0.0) # >0 means filter by NaNs
    
    # Read full expression
    df = pd.read_csv(expression_file, sep='\t', index_col=0)
    
    # Subset samples
    # Ensure intersection first
    common = [s for s in samples if s in df.columns]
    if len(common) < len(samples):
        print(f"Warning: Only {len(common)}/{len(samples)} requested samples found in expression matrix.")
    df = df[common]
    
    # Filter
    if filter_low_expr:
        initial_count = len(df)
        
        # 1. Filter by Low Expression
        skip_tpm = filter_args.get('skip_tpm', False)
        
        n_samples = len(common)
        
        if not skip_tpm:
            # Calc threshold count
            thresh = min_samples_abs if min_samples_abs > 0 else int(n_samples * min_samples_rel)
            
            # Keep if > min_tpm in >= thresh samples
            mask_expr = (df > min_tpm).sum(axis=1) >= thresh
            df = df[mask_expr]
            print(f"Expression Filter (>{min_tpm} in >={thresh} samples): {initial_count} -> {len(df)} genes.")
        else:
            print("Skipping low expression value filter (User requested).")
        
        # 2. Filter by Missing (NaN)
        if max_missing > 0:
            # Drop if nan_count / n_samples > max_missing
            current_count = len(df)
            nan_frac = df.isna().sum(axis=1) / n_samples
            df = df[nan_frac <= max_missing]
            print(f"Missing Filter (NaN <= {max_missing*100}%): {current_count} -> {len(df)} genes.")
            
        # Fill NA with 0 before INT? Or handle in INT?
        # Usually for INT we impute or drop. Let's fill 0 for now as safely assumed zero expr.
        df = df.fillna(0)
        
    # Write temp for IQN
    temp_file = output_dir / "expression_filtered.tsv"
    df.to_csv(temp_file, sep='\t')
    
    final_file = temp_file
    
    # Normalization (IQN) using existing script logic
    if normalize_iqn:
        iqn_file = output_dir / "expression.qnorm"
        # We can implement IQN in python or call the R script
        # For simplicity/speed, implementing basic Rank based INT in python
        # or call the moved script `analysis/run_iqn.py`
        
        # Python implementation of INT (Inverse Normal Transformation)
        from scipy.stats import rankdata, norm
        
        def inverse_normal_transform(series):
            # Rank, average ties
            ranks = rankdata(series, method='average')
            # Convert to quantiles
            quantiles = (ranks - 0.5) / len(ranks)
            # Inverse CDF
            return norm.ppf(quantiles)

        # Apply row-wise (per gene)
        # Apply row-wise (per gene)
        # Using result_type='expand' to ensure we get a DataFrame back, not a Series of arrays
        df_norm = df.apply(inverse_normal_transform, axis=1, result_type='expand')
        
        # Ensure columns are preserved (though apply(axis=1) usually preserves index as index, and result becomes columns)
        df_norm.columns = df.columns
        df_norm.to_csv(iqn_file, sep='\t')
        final_file = iqn_file
        
    return str(final_file)

    return str(final_file)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Expression Preprocessing (Filter + IQN)")
    parser.add_argument("--expression", required=True, help="Input expression matrix file")
    parser.add_argument("--samples", required=True, help="File with list of samples to keep")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    
    # Filter args
    parser.add_argument("--no-filter", action="store_true", help="Skip low expression filtering")
    parser.add_argument("--no-tpm-filter", action="store_true", help="Skip value-based filtering (TPM/Threshold), only do Missing filter if requested")
    parser.add_argument("--min-tpm", type=float, default=0.1, help="Minimum expression value (TPM) to consider 'detected'")
    parser.add_argument("--min-samples", type=int, default=0, help="Minimum number of samples with detected expression")
    parser.add_argument("--min-samples-rel", type=float, default=0.2, help="Minimum proportion of samples (0.0-1.0) if --min-samples not set")
    parser.add_argument("--max-missing", type=float, default=0.0, help="Filter genes with > proportion of missing (NaN) values")
    
    parser.add_argument("--no-norm", action="store_true", help="Skip IQN normalization")
    
    args = parser.parse_args()
    
    # Read samples
    with open(args.samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]

    # Add flexible logic to process_expression (need to update signature too, but for CLI we can pass kwargs or just modify function above)
    # Let's update process_expression signature first.
        
    process_expression(
        expression_file=args.expression, 
        samples=samples, 
        output_dir=args.out_dir,
        filter_low_expr=not args.no_filter,
        normalize_iqn=not args.no_norm,
        filter_args={
            'min_tpm': args.min_tpm,
            'min_samples': args.min_samples,
            'min_samples_rel': args.min_samples_rel,
            'max_missing': args.max_missing,
            'skip_tpm': args.no_tpm_filter
        }
    )

if __name__ == "__main__":
    main()
