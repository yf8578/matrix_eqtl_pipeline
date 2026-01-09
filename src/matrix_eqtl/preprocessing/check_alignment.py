import sys
import os

def get_samples(filepath):
    """Read the header of a file and return sample IDs (skipping first column)."""
    if not os.path.exists(filepath):
        print(f"Error: File not found: {filepath}")
        sys.exit(1)
    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        # Assuming first column is ID (SNP/Gene/Cov), rest are samples
        if len(header) < 2:
            return []
        return header[1:]

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 check_alignment.py [GenotypeFile] [ExpressionFile] [CovariatesFile]")
        sys.exit(1)

    geno_file = sys.argv[1]
    expr_file = sys.argv[2]
    cov_file = sys.argv[3]

    print("Checking Sample Alignment...")
    print(f"1. Genotype:   {os.path.basename(geno_file)}")
    print(f"2. Expression: {os.path.basename(expr_file)}")
    print(f"3. Covariates: {os.path.basename(cov_file)}")
    print("-" * 40)

    s1 = get_samples(geno_file)
    s2 = get_samples(expr_file)
    s3 = get_samples(cov_file)

    print(f"Sample Counts -> G: {len(s1)}, E: {len(s2)}, C: {len(s3)}")

    # Check 1: Counts
    if not (len(s1) == len(s2) == len(s3)):
        print("❌ FAIL: Sample counts do not match!")
        sys.exit(1)

    # Check 2: Exact Match
    mismatch = False
    if s1 != s2:
        print("❌ FAIL: Genotype and Expression samples differ!")
        print(f"   First 5 G: {s1[:5]}")
        print(f"   First 5 E: {s2[:5]}")
        mismatch = True
    
    if s1 != s3:
        print("❌ FAIL: Genotype and Covariate samples differ!")
        print(f"   First 5 G: {s1[:5]}")
        print(f"   First 5 C: {s3[:5]}")
        mismatch = True

    if not mismatch:
        print("✅ SUCCESS: All samples are perfectly aligned!")

if __name__ == "__main__":
    main()
