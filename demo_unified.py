
import os
import sys
from pathlib import Path
from matrix_eqtl.pipeline import EQTLPipeline

# Setup paths (Assuming running from project root)
PROJECT_ROOT = Path(".").resolve()
EXAMPLES_DIR = PROJECT_ROOT / "examples"

# These files were previously created by verify_installation or create_demo_data
# We can re-use them if they exist, or create new simple ones.
# Let's assume we run create_demo_data.py first.

def create_dummy_data():
    import numpy as np
    import pandas as pd
    
    data_dir = EXAMPLES_DIR / "input"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    samples = [f"Sample{i}" for i in range(1, 21)]
    
    # Expression
    expr = pd.DataFrame(np.random.rand(100, 20), index=[f"Gene{i}" for i in range(100)], columns=samples)
    expr_file = data_dir / "expression.tsv"
    expr.to_csv(expr_file, sep='\t')
    
    # VCF
    vcf_file = data_dir / "genotypes.vcf"
    with open(vcf_file, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for i in range(1, 101):
            # Simple bi-allelic SNPs
            gt = ["0/0", "0/1", "1/1"]
            gts = [gt[np.random.randint(0,3)] for _ in range(20)]
            f.write(f"1\t{i*1000}\trs{i}\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
            
    # Covariates
    cov = pd.DataFrame(np.random.rand(2, 20), index=["Age", "Sex"], columns=samples)
    cov_file = data_dir / "covariates.txt"
    cov.to_csv(cov_file, sep='\t')
    
    # Gene/SNP Loc
    gene_loc = pd.DataFrame({
        'geneid': [f"Gene{i}" for i in range(100)],
        'chr': ['1'] * 100,
        'left': [i*1000 for i in range(1, 101)],
        'right': [i*1000 + 100 for i in range(1, 101)]
    })
    gene_loc_file = data_dir / "gene_loc.txt"
    gene_loc.to_csv(gene_loc_file, sep='\t', index=False)
    
    return str(vcf_file), str(expr_file), str(cov_file), str(gene_loc_file)

def main():
    print("Generating dummy data...")
    vcf, expr, cov, gene_loc = create_dummy_data()
    
    # Initialize Pipeline
    output_dir = PROJECT_ROOT / "examples" / "output_demo"
    
    # Test MatrixEQTL
    print("\n>>> Testing MatrixEQTL...")
    pipeline = EQTLPipeline(output_dir=output_dir / "matrix_eqtl", tool='matrixeqtl')
    pipeline.run(
        vcf_file=vcf,
        expression_file=expr,
        covariates_file=cov,
        gene_location_file=gene_loc,
        genotype_pca_n=2,
        peer_factors_n=0, # Need PEER installed for this to work, skip for basic test
        cis_window=5000
    )
    
    # Test QTLtools
    # Note: QTLtools binary must be in PATH. If not, this will fail.
    # We can skip if check fails.
    if os.system("which QTLtools > /dev/null") == 0:
        print("\n>>> Testing QTLtools...")
        pipeline_q = EQTLPipeline(output_dir=output_dir / "qtltools", tool='qtltools')
        pipeline_q.run(
            vcf_file=vcf,
            expression_file=expr,
            covariates_file=cov,
            gene_location_file=gene_loc,
            genotype_pca_n=2,
            cis_window=5000
        )
    else:
        print("\n>>> QTLtools not found, skipping test.")

if __name__ == "__main__":
    main()
