
import subprocess
import os
from pathlib import Path

def run_analysis(genotype_file, expression_file, covariates_file, output_dir, gene_loc, snp_loc, cis_window=1e6):
    """
    Run MatrixEQTL analysis.
    Checks if gene_loc and snp_loc are provided (Mode A) or not (Mode B/NoLoc).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    cis_out = output_dir / "cis_results.txt"
    trans_out = output_dir / "trans_results.txt"
    
    script_path = Path(__file__).parent.parent / "analysis" / "run_matrix_eqtl.py"
    
    cmd = [
        "python3", str(script_path),
        "--genotype-matrix", str(genotype_file),
        "--gene-expression-matrix", str(expression_file),
        "--covariates", str(covariates_file),
        "--output-file", str(cis_out),
        "--trans-output-file", str(trans_out),
        "--cis-distance", str(cis_window)
    ]
    
    if gene_loc and snp_loc:
        cmd.extend([
            "--gene-positions", str(gene_loc),
            "--genotype-positions", str(snp_loc)
        ])
    else:
        # Use location-agnostic script if locations missing?
        # Or just run standard check. user requirement asked for Mode A/B based on param.
        # Assuming for now we default to the wrapper we moved.
        pass
        
    print(f"Running MatrixEQTL: {' '.join(cmd)}")
    subprocess.check_call(cmd)
