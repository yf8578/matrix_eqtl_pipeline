import subprocess
import argparse
import sys
import os
import pandas as pd
from check_file import check_file

def get_args():
    parser = argparse.ArgumentParser(description="Process VCF using PLINK for MatrixEQTL")
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF file (supports .vcf.gz)')
    parser.add_argument('-e', '--expression', required=True, help='Expression file to extract sample IDs')
    parser.add_argument('-o', '--out-dir', required=True, help='Output directory')
    parser.add_argument('--maf', type=float, default=0.05, help='MAF threshold (default: 0.05)')
    parser.add_argument('--plink-cmd', default='plink', help='Command to invoke PLINK (default: plink)')
    
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    
    # Ensure output directory exists
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    print(f"Step 1: Extracting sample IDs from {args.expression}...")
    # Read expression file header (tab-separated)
    try:
        with open(args.expression, 'r') as f:
            header = f.readline().strip().split('\t')
            # Assuming first column is GeneID, rest are Samples
            samples = header[1:]
            
        print(f"  Found {len(samples)} samples.")
        
        # Write format for PLINK --keep (FamilyID SampleID)
        # Since VCF usually has FID=SampleID or just SampleID
        # We assume FID = 0 or FID = IID depending on VCF. 
        # For standard VCF imported to PLINK, usually --double-id or --const-fid is handled,
        # but --keep accepts IID lists if we format correctly.
        # Safe bet: Write SampleID SampleID (if PLINK uses standard import where FID=IID)
        # Or just write SampleID if we use --keep <filename> with one column? 
        # PLINK 1.9 --keep usually expects 2 cols: FID IID.
        # Let's try to be robust. Most standard VCF tools map Sample ID to both FID and IID.
        keep_file = os.path.join(args.out_dir, "samples_to_keep.txt")
        with open(keep_file, 'w') as f:
            for s in samples:
                f.write(f"{s}\t{s}\n")
                
    except Exception as e:
        print(f"Error reading expression file: {e}")
        sys.exit(1)

    print(f"Step 2: Running PLINK...")
    plink_prefix = os.path.join(args.out_dir, "plink_temp")
    
    # Construct PLINK command
    # --vcf: input
    # --keep: subset samples
    # --maf: filter variants
    # --recode A-transpose: output .traw format (raw counts, optimal for MatrixEQTL)
    # --allow-no-sex: prevent errors on missing sex info
    # --double-id: set FID = IID (crucial for VCFs without family info)
    cmd = [
        args.plink_cmd,
        "--vcf", args.vcf,
        "--keep", keep_file,
        "--maf", str(args.maf),
        "--recode", "A-transpose",
        "--allow-no-sex",
        "--double-id",
        "--out", plink_prefix
    ]
    
    print("  Command: " + " ".join(cmd))
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        print("PLINK failed. Please check if 'plink' is installed and accessible.")
        sys.exit(1)
    
    print(f"Step 3: Formatting output for MatrixEQTL...")
    traw_file = plink_prefix + ".traw"
    final_matrix = os.path.join(args.out_dir, "genotype.matrix")
    final_pos = os.path.join(args.out_dir, "snp_positions.txt")
    
    # .traw format: CHR SNP  (cM) POS COUNTED ALLELE1 ALLELE2 SAMPLE1 SAMPLE2 ...
    # We need:
    # 1. Genotype Matrix: SNPID Sample1 Sample2 ...
    # 2. Position File: SNPID chr pos
    
    # We use pandas for efficient processing (chunking if really huge, but read_csv is usually fine for <100k SNPs)
    # If huge, use awk. Let's start with pandas for code simplicity as requested.
    
    try:
        # Read .traw
        df = pd.read_csv(traw_file, sep='\t')
        
        # 3.1 Positions: SNP, CHR, POS
        # Ensure CHR is just numeric (PLINK usually handles this, ensures 1 not chr1)
        pos_df = df[['SNP', 'CHR', 'POS']].copy()
        pos_df.columns = ['snpid', 'chr', 'pos']
        pos_df.to_csv(final_pos, sep='\t', index=False)
        
        # 3.2 Matrix: ID, Samples...
        # Drop metadata cols: CHR, (cM), POS, COUNTED, ALLELE1, ALLELE2
        # Traw cols: CHR, SNP, (cM), POS, COUNTED, ALT, ...Samples...
        # Note: PLINK 1.9 output usually has CHR SNP (C)M POS COUNTED ALT ...
        # Let's drop by name to be safe
        cols_to_drop = ['CHR', 'POS', 'COUNTED', 'ALT', '(C)M']
        # If (C)M not present (version diff), ignore errors
        matrix_df = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
        
        # Rename SNP to id (or whatever MatrixEQTL expects, usually just ID in first col)
        matrix_df.rename(columns={'SNP': 'id'}, inplace=True)
        
        # Clean Sample IDs (PLINK .traw format: FID_IID or just IID depending on loading)
        # Since we used --double-id, they are likely Sample_Sample.
        # But wait, .traw usually just puts IID if unique?
        # Let's clean the columns: split by '_' and take last part IF they look like dups
        # Actually simplest verification: Check with user sample list
        current_cols = matrix_df.columns.tolist()
        new_cols = []
        for c in current_cols:
            if c == 'id':
                new_cols.append(c)
            else:
                # If format is FID_IID and FID=IID, it looks like SAMPLE_SAMPLE
                if f"{c}_{c}" in current_cols: 
                    # This logic is tricky. Let's just strip the prefix if it matches
                    # format: ID_ID
                    if '_' in c:
                        parts = c.split('_')
                        if len(parts) == 2 and parts[0] == parts[1]:
                             new_cols.append(parts[0])
                        else:
                             new_cols.append(c)
                    else:
                        new_cols.append(c)
        matrix_df.columns = new_cols
        
        # Save
        matrix_df.to_csv(final_matrix, sep='\t', index=False)
        
        print(f"  Genotype Matrix: {final_matrix}")
        print(f"  SNP Positions:   {final_pos}")

        # Cleanup
        os.remove(keep_file)
        # Maybe keep PLINK logs? User asked "will I use decompressed VCF", we say no.
        # We can rm the .traw to save space
        os.remove(traw_file)

    except Exception as e:
        print(f"Error processing PLINK output: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
