import argparse
import subprocess
import os
import sys
import pandas as pd
def get_args():
    parser = argparse.ArgumentParser(description='Process VCF using PLINK for MatrixEQTL')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--expression', required=True, help='Expression matrix (to extract sample IDs)')
    parser.add_argument('-o', '--out-dir', required=True, help='Output directory')
    parser.add_argument('--maf', type=float, default=0.05, help='MAF threshold (default: 0.05)')
    parser.add_argument('--plink-cmd', default='plink', help='Command to invoke PLINK (default: plink)')
    parser.add_argument('--snps-only', action='store_true', help='Only analyze SNPs (exclude indels)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for PLINK')
    
    args = parser.parse_args()
    return args
def main():
    args = get_args()
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    print(f"Step 1: Extracting sample IDs from {args.expression}...")
    try:
        with open(args.expression, 'r') as f:
            header = f.readline().strip().split('\t')
            samples = header[1:] # Assuming first col is GeneID
        print(f"  Found {len(samples)} samples.")
        
        keep_file = os.path.join(args.out_dir, "samples_to_keep.txt")
        with open(keep_file, 'w') as f:
            for s in samples:
                f.write(f"{s}\t{s}\n")
    except Exception as e:
        print(f"Error reading expression file: {e}")
        sys.exit(1)
    print(f"Step 2: Running PLINK (Threads: {args.threads})...")
    plink_prefix = os.path.join(args.out_dir, "plink_temp")
    
    cmd = [
        args.plink_cmd,
        "--vcf", args.vcf,
        "--keep", keep_file,
        "--maf", str(args.maf),
        "--make-bed",          
        "--recode", "A-transpose", 
        "--allow-no-sex",
        "--double-id",
        "--out", plink_prefix,
        "--threads", str(args.threads)
    ]
    if args.snps_only:
        cmd.append("--snps-only")
    
    # print("  Command: " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"PLINK failed with error code {e.returncode}")
        sys.exit(1)
    print(f"Step 3: Formatting output for MatrixEQTL (Optimized Chunked Processing)...")
    
    traw_file = plink_prefix + ".traw"
    final_matrix = os.path.join(args.out_dir, "genotype.matrix")
    final_pos = os.path.join(args.out_dir, "snp_positions.txt")

    # 1. Process SNP Positions (Fast)
    # Read specific columns for positions
    print("  Processing positions...")
    try:
        # Usecols: 0=CHR, 1=SNP, 3=POS. (cM is col 2, ignored)
        pos_df = pd.read_csv(traw_file, sep='\t', usecols=['CHR', 'SNP', 'POS'])
        # Rename to MatrixEQTL format: snpid, chr, pos
        pos_df = pos_df[['SNP', 'CHR', 'POS']]
        pos_df.columns = ['snpid', 'chr', 'pos']
        # Add 'chr' prefix if needed, usually PLINK output is just number
        # MatrixEQTL handles both usually, but let's keep it as is.
        pos_df.to_csv(final_pos, sep='\t', index=False)
        print(f"  Saved positions to {final_pos}")
    except Exception as e:
        print(f"Error processing positions: {e}")
        sys.exit(1)

    # 2. Process Genotype Matrix (Chunked)
    print("  Processing genotype matrix (Chunked)...")
    chunk_size = 100000 
    first_chunk = True
    
    try:
        # We need to skip the first 6 columns metadata to get just the genotype data
        # But we need row ID from 'SNP' column.
        # So we read everything but process in chunks.
        
        # Read iterator
        reader = pd.read_csv(traw_file, sep='\t', chunksize=chunk_size)
        
        with open(final_matrix, 'w') as f_out:
            for chunk in reader:
                # Store SNP column as index for this chunk
                chunk.index = chunk['SNP']
                
                # Drop metadata columns: CHR, (SNP), (cM), POS, COUNTED, ALT
                # Traw standard: CHR, SNP, (cM), POS, COUNTED, ALT, ...Samples...
                # So we drop first 6 columns.
                data_chunk = chunk.iloc[:, 6:]
                
                # Fix Column headers (Sample IDs)
                if first_chunk:
                    # Clean headers: remove FID_ prefix if present (0_Sample1 -> Sample1)
                    new_cols = []
                    for c in data_chunk.columns:
                        if '_' in c:
                            parts = c.split('_', 1)
                            # Simple heuristic: if parts[0] is '0', take parts[1]
                            # Or if parts[0] == parts[1] (double-id), take parts[0]
                            if parts[0] == parts[1]: 
                                new_cols.append(parts[0])
                            elif parts[0] == '0':
                                new_cols.append(parts[1])
                            else:
                                new_cols.append(c) # Fallback
                        else:
                            new_cols.append(c)
                    data_chunk.columns = new_cols
                
                # Write to file
                # mode='a' if not first, but we opened file 'w' outside loop so we just write
                # header only for first chunk
                data_chunk.to_csv(f_out, sep='\t', index=True, index_label='id', header=first_chunk, na_rep='NA', float_format='%.0f')
                
                first_chunk = False
                print(f"    Processed {chunk.shape[0]} rows...", end='\r')
                
        print(f"\n  Saved genotype matrix to {final_matrix}")
        
    except Exception as e:
        print(f"Error processing matrix: {e}")
        sys.exit(1)
if __name__ == '__main__':
    main()