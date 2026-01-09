
import pandas as pd
import gzip

def get_common_samples(vcf_file: str, expression_file: str) -> list:
    """
    Find common samples between VCF and Expression matrix.
    """
    # 1. Get VCF samples
    vcf_samples = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                vcf_samples = line.strip().split('\t')[9:]
                break
                
    # 2. Get Expression samples
    # Assuming expression file has samples as columns
    df = pd.read_csv(expression_file, sep='\t', index_col=0, nrows=1)
    exp_samples = df.columns.tolist()
    
    # 3. Intersect
    common = sorted(list(set(vcf_samples) & set(exp_samples)))

    if not common:
        print("ERROR: No common samples found between VCF and Expression file.")
        print(f"  VCF Total: {len(vcf_samples)} samples. First 5: {vcf_samples[:5]}")
        print(f"  Expression Total: {len(exp_samples)} samples. First 5: {exp_samples[:5]}")
        raise ValueError("Sample intersection is empty. Check your sample IDs.")
        
    return common

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Sample Intersection")
    parser.add_argument("--vcf", required=True, help="VCF file")
    parser.add_argument("--expression", required=True, help="Expression file")
    parser.add_argument("--out", required=True, help="Output file for common samples list")
    
    args = parser.parse_args()
    
    samples = get_common_samples(args.vcf, args.expression)
    print(f"Found {len(samples)} common samples.")
    
    with open(args.out, 'w') as f:
        for s in samples:
            f.write(s + "\n")

if __name__ == "__main__":
    main()
