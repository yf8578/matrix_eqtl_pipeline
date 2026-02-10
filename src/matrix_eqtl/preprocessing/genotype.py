
import subprocess
from pathlib import Path
import os

def process_genotype(vcf_file, samples, output_dir, target_tool, maf="0.01", geno="0.05", hwe="1e-6", threads=1, plink_bin="plink2", plink_version=2):
    """
    Process VCF using PLINK to generate required format.
    plink_version: 1 (for v1.9) or 2 (for v2.0)
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
    
    # PLINK output prefix
    plink_out = output_dir / "genotypes"
    
    # Run PLINK QC
    if plink_version == 2:
        # PLINK 2 cannot --make-bed --sort-vars directly.
        # It requires: VCF -> Sorted PGEN -> BED.
        
        # Step 1: VCF -> Sorted PGEN
        temp_pgen_prefix = output_dir / "temp_sorted"
        cmd_step1 = [
            plink_bin,
            "--vcf", str(vcf_file),
            "--keep", str(keep_file),
            "--make-pgen",
            "--sort-vars",
            "--out", str(temp_pgen_prefix),
            "--maf", str(maf),
            "--geno", str(geno),
            "--hwe", str(hwe),
            "--threads", str(threads),
            "--double-id",
            "--allow-extra-chr"
        ]
        print(f"Running PLINK 2 Step 1 (Sort & QC): {' '.join(cmd_step1)}")
        subprocess.check_call(cmd_step1)
        
        # Step 2: Sorted PGEN -> BED (Fixed-width)
        cmd_step2 = [
            plink_bin,
            "--pfile", str(temp_pgen_prefix),
            "--make-bed",
            "--out", str(plink_out),
            "--threads", str(threads)
        ]
        print(f"Running PLINK 2 Step 2 (PGEN -> BED): {' '.join(cmd_step2)}")
        subprocess.check_call(cmd_step2)
        
        # Clean up temp PGEN files
        for ext in ['.pgen', '.pvar', '.psam', '.log']:
            f_rem = Path(str(temp_pgen_prefix) + ext)
            if f_rem.exists():
                f_rem.unlink()
                
    else:
        # PLINK 1.9 Workflow (Direct)
        cmd = [
            plink_bin,
            "--vcf", str(vcf_file),
            "--keep", str(keep_file),
            "--make-bed",
            "--out", str(plink_out),
            "--maf", str(maf),
            "--geno", str(geno),
            "--hwe", str(hwe),
            "--threads", str(threads),
            "--double-id",
            "--allow-extra-chr",
            "--allow-no-sex"
        ]
        print(f"Running PLINK QC (v{plink_version}): {' '.join(cmd)}")
        subprocess.check_call(cmd)
    
    result = {
        'plink_prefix': str(plink_out)
    }
    
    if target_tool == 'matrixeqtl':
        # Use --export A-transpose (PLINK 2 syntax)
        # Output format is .traw
        export_prefix = output_dir / "genotypes_export"
        cmd_recode = []
        
        if plink_version == 2:
            # Use --export A-transpose (PLINK 2 syntax)
            cmd_recode = [
                plink_bin, "--bfile", str(plink_out), 
                "--export", "A-transpose", 
                "--out", str(export_prefix)
            ]
        else:
            # Use --recode A-transpose (PLINK 1.9 syntax)
            cmd_recode = [
                plink_bin, "--bfile", str(plink_out), 
                "--recode", "A-transpose", 
                "--out", str(export_prefix)
            ]
            
        print(f"Running PLINK 2 Export (Transpose): {' '.join(cmd_recode)}")
        subprocess.check_call(cmd_recode)
        
        traw_file = str(export_prefix) + ".traw"
        matrix_file = output_dir / "genotype.matrix"
        
        print("Formatting genotype matrix (chunked processing)...")
        
        # .traw headers: CHR SNP (CM) POS COUNTED ALLELE FID_IID ...
        # matrixeqtl wants: id sample1 sample2 ...
        
        chunk_size = 50000 
        # Metadata cols to drop (PLINK .traw format variances)
        drop_cols = ['CHR', '(CM)', '(C)M', 'POS', 'COUNTED', 'ALLELE', 'ALT']
        
        first_chunk = True
        
        import pandas as pd
        for chunk in pd.read_csv(traw_file, sep='\t', chunksize=chunk_size, dtype={'CHR': str}):
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
                # ... (rest of logic remains same but truncated in replacement for brevity if needed, but here I must provide full block or matching end)
                for c in chunk.columns:
                     # ... (re-implementing the inner logic to match exactly what is there or use context)
                     if c == 'id': continue
                     
                     clean_id = c
                     if c not in samples:
                          # Try removing duplication "X_X" -> "X" caused by FID=IID
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
        bim = pd.read_csv(bim_file, sep='\t', header=None, names=['chr','id','cm','pos','a1','a2'], dtype={'chr': str})
        bim[['id','chr','pos']].to_csv(snp_loc_file, sep='\t', index=False)
        result['snp_loc'] = str(snp_loc_file)

    elif target_tool == 'qtltools':
        # QTLtools needs VCF.
        # Use PLINK 2 to export VCF with filtered samples/SNPs
        vcf_out = output_dir / "genotypes.filtered.vcf.gz"
        cmd_vcf = []
        
        if plink_version == 2:
            # Use --export vcf bgz (PLINK 2 syntax)
            # Add id-paste=iid to prevent PLINK from munging IDs (e.g. FID_IID)
            # This ensures VCF IDs match the input sample list exactly
            cmd_vcf = [
                plink_bin, "--bfile", str(plink_out), 
                "--export", "vcf", "bgz", "id-paste=iid", 
                "--out", str(output_dir / "genotypes.filtered")
            ]
        else:
            # Use --recode vcf-iid bgz (PLINK 1.9 syntax)
            cmd_vcf = [
                plink_bin, "--bfile", str(plink_out), 
                "--recode", "vcf-iid", "bgz", 
                "--out", str(output_dir / "genotypes.filtered")
            ]
            
        subprocess.check_call(cmd_vcf)
        
        # Output handling
        # plink2 --export vcf bgz -> .vcf.gz
        # plink1 --recode vcf-iid bgz -> .vcf.gz (usually)
        # Check generated file
        generated = output_dir / "genotypes.filtered.vcf.gz"
        if not generated.exists():
             # Try without gz extension if plink1 behavior differs or find what it made
             pass

        # Rename output
        os.rename(str(generated), str(vcf_out))
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
    
    # PLINK Version Control
    parser.add_argument("--plink-bin", default="plink2", help="Path or name of PLINK executable (default: plink2)")
    parser.add_argument("--plink-version", type=int, default=2, choices=[1, 2], help="PLINK version: 1 for v1.9, 2 for v2.0 (default: 2)")
    
    args = parser.parse_args()
    
    # Read samples
    with open(args.samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
        
    process_genotype(args.vcf, samples, args.out_dir, args.tool, 
                     maf=args.maf, geno=args.geno, hwe=args.hwe, threads=args.threads,
                     plink_bin=args.plink_bin, plink_version=args.plink_version)

if __name__ == "__main__":
    main()