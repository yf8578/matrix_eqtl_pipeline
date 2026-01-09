
import subprocess
from pathlib import Path
import os

def process_genotype(vcf_file, samples, output_dir, target_tool, maf="0.01", geno="0.05", hwe="1e-6", threads=1):
    """
    Process VCF using PLINK to generate required format.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write sample list
    keep_file = output_dir / "keep_samples.txt"
    with open(keep_file, 'w') as f:
        for s in samples:
            f.write(f"{s} {s}\n") # PLINK expects FID IID
            
    # PLINK output prefix
    plink_out = output_dir / "genotypes"
    
    # Run PLINK: Filter samples, make bed
    cmd = [
        "plink",
        "--vcf", str(vcf_file),
        "--keep", str(keep_file),
        "--make-bed",
        "--out", str(plink_out),
        "--maf", str(maf),
        "--geno", str(geno),
        "--hwe", str(hwe),
        "--threads", str(threads)
    ]
    # Add --const-fid if VCF only has IID? Assuming standard VCF.
    # Often VCF sample IDs are treated as IID, FID=0 or FID=IID. 
    # Let's assume standard behavior first.
    
    print(f"Running PLINK QC: {' '.join(cmd)}")
    subprocess.check_call(cmd)
    
    result = {
        'plink_prefix': str(plink_out)
    }
    
    if target_tool == 'matrixeqtl':
        # Use --recode A-transpose to get SNPs in rows (Variant-major)
        # This avoids loading a massive (Samples x SNPs) matrix into pandas to transpose it.
        # Output format is .traw
        # Use a distinct prefix for export to avoid any potential conflict/overwrite of the main BED
        export_prefix = output_dir / "genotypes_export"
        cmd_recode = [
            "plink", "--bfile", str(plink_out), 
            "--recode", "A-transpose", 
            "--out", str(export_prefix)
        ]
        print(f"Running PLINK Recode (Transpose): {' '.join(cmd_recode)}")
        subprocess.check_call(cmd_recode)
        
        traw_file = str(export_prefix) + ".traw"
        matrix_file = output_dir / "genotype.matrix"
        
        print("Formatting genotype matrix (chunked processing)...")
        
        # .traw headers: CHR SNP (CM) POS COUNTED ALLELE FID_IID ...
        # matrixeqtl wants: id sample1 sample2 ...
        
        chunk_size = 50000 
        # Metadata cols to drop
        drop_cols = ['CHR', '(CM)', 'POS', 'COUNTED', 'ALLELE']
        
        first_chunk = True
        
        import pandas as pd
        for chunk in pd.read_csv(traw_file, sep='\t', chunksize=chunk_size):
            # 1. Drop metadata
            chunk.drop(columns=[c for c in drop_cols if c in chunk.columns], inplace=True)
            
            # 2. Rename 'SNP' -> 'id' (optional, MatrixEQTL handles first col as ID usually)
            if 'SNP' in chunk.columns:
                chunk.rename(columns={'SNP': 'id'}, inplace=True)
                
            # 3. Clean Sample IDs (FID_IID -> IID)
            # Since we used "s s" for keep file, PLINK usually outputs "s_s"
            # We need to map them back to "s"
            # We can detect this from the header of the first chunk
            if first_chunk:
                new_cols = {}
                for c in chunk.columns:
                    if c == 'id': continue
                    # Logic: if column looks like "SampleA_SampleA", take "SampleA"
                    # Or just take the part after the last underscore if we are sure?
                    # Safer: known samples check?
                    # Let's assume FID_IID format.
                    # Since we set FID=IID, we can just split by _ and take the last part? 
                    # Or better: check if `s_s` is in the samples list we want? No, we want `s`.
                    # Let's define a cleaner function:
                    if '_' in c:
                        # try exact match with our samples list?
                        # This might be slow for 10M rows loop, but we only do it on columns map
                        pass
                    
                    # Heuristic: split on first `_`?
                    # If ID itself has `_`, PLINK joins FID_IID with `_`.
                    # If FID=IID="A_B", PLINK might output "A_B_A_B".
                    # If we set FID=IID, we expect "ID_ID".
                    # Let's try to match the prefix.
                    
                    # Robust cleanup: 
                    # If col starts with its own suffix + "_"?
                    # Actually, we passed `samples`. We know what we entered.
                    # We just need to find which col matches which sample.
                    # PLINK preserves simple IDs usually.
                    
                    # Simple fix: split by LAST underscore if FID=IID
                    # But if ID has no underscore, PLINK might output 0_ID?
                    # We used `keep_samples.txt` with `s s`. So FID=s.
                    # So PLINK outputs `s_s`.
                    
                    # Strategy: If ID is in our expected samples, keep it.
                    # If not, try stripping `s_`.
                    # Pythonic:
                    clean_id = c
                    if c not in samples:
                         # Try removing duplication "X_X" -> "X"
                         # Check if string is two halves equal?
                         mid = len(c) // 2
                         if len(c) % 2 == 1 and c[mid] == '_':
                             if c[:mid] == c[mid+1:]:
                                 clean_id = c[:mid]
                    
                    new_cols[c] = clean_id
                    
                chunk.rename(columns=new_cols, inplace=True)
                
                # Verify intersection with requested samples
                # (Optional warning)
                
                chunk.to_csv(matrix_file, sep='\t', index=False, mode='w')
                first_chunk = False
            else:
                chunk.to_csv(matrix_file, sep='\t', index=False, mode='a', header=False)

        result['matrix'] = str(matrix_file)
        
        # SNP Loc
        # .bim file is still generated by --make-bed earlier
        bim_file = str(plink_out) + ".bim"
        snp_loc_file = output_dir / "snp_positions.txt"
        bim = pd.read_csv(bim_file, sep='\t', header=None, names=['chr','id','cm','pos','a1','a2'])
        bim[['id','chr','pos']].to_csv(snp_loc_file, sep='\t', index=False)
        result['snp_loc'] = str(snp_loc_file)

    elif target_tool == 'qtltools':
        # QTLtools needs VCF.
        # It's already VCF, but might need to be subsetted/re-ordered.
        # Use PLINK to export VCF with filtered samples/SNPs
        vcf_out = output_dir / "genotypes.filtered.vcf.gz"
        cmd_vcf = [
            "plink", "--bfile", str(plink_out), 
            "--recode", "vcf-iid", "bgz", 
            "--out", str(output_dir / "genotypes.filtered")
        ]
        subprocess.check_call(cmd_vcf)
        
        # Rename output
        os.rename(str(output_dir / "genotypes.filtered.vcf.gz"), str(vcf_out))
        # Index it
        subprocess.check_call(["tabix", "-p", "vcf", str(vcf_out)])
        
        result['vcf'] = str(vcf_out)
        
    return result

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Genotype Preprocessing (PLINK wrapper)")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--samples", required=True, help="File with list of samples to keep (one per line)")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--tool", required=True, choices=['matrixeqtl', 'qtltools', 'tensorqtl'], help="Target analysis tool")
    
    # PLINK Params
    parser.add_argument("--maf", default="0.01", help="MAF threshold (default 0.01)")
    parser.add_argument("--geno", default="0.05", help="Geno missing threshold (default 0.05)")
    parser.add_argument("--hwe", default="1e-6", help="HWE threshold (default 1e-6)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for PLINK")
    # Optional strict args?
    
    args = parser.parse_args()
    
    # Read samples
    with open(args.samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
        
    process_genotype(args.vcf, samples, args.out_dir, args.tool, 
                     maf=args.maf, geno=args.geno, hwe=args.hwe, threads=args.threads)

if __name__ == "__main__":
    main()