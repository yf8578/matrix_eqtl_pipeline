
import subprocess
from pathlib import Path

def run_analysis(vcf_file, expression_bed, covariates_file, output_dir, gene_loc, cis_window=1000000):
    """
    Run QTLtools cis analysis.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "qtltools_cis.txt"
    
    # Check for qtltools binary (assuming user has it in path or conda)
    qtltools_cmd = "QTLtools" 
    
    cmd = [
        qtltools_cmd, "cis",
        "--vcf", str(vcf_file),
        "--bed", str(expression_bed),
        "--cov", str(covariates_file),
        "--out", str(output_file),
        "--window", str(cis_window),
        "--nominal", "1" # or --permute 1000
    ]
    
    print(f"Running QTLtools: {' '.join(cmd)}")
    subprocess.check_call(cmd)
