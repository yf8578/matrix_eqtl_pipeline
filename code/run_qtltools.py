import argparse
import os
import subprocess
import sys

def check_qtltools():
    """Check if QTLtools is in PATH."""
    try:
        subprocess.check_call(['which', 'QTLtools'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("Error: QTLtools not found in PATH.")
        print("Please install it or export PATH=/opt/qtltools/bin:$PATH")
        sys.exit(1)

def get_args():
    parser = argparse.ArgumentParser(description="Wrapper for QTLtools Analysis. \nAny extra arguments (e.g. --seed 123) provided after the script arguments will be passed directly to QTLtools.")
    
    # Input Files
    parser.add_argument('--vcf', required=True, help='Genotype VCF/BCF file')
    parser.add_argument('--bed', required=True, help='Phenotype BED file (bgzipped & indexed)')
    parser.add_argument('--cov', help='Covariates file (optional)')
    
    # Output
    parser.add_argument('--out', required=True, help='Output file path')
    
    # Mode
    parser.add_argument('--mode', choices=['permute', 'nominal', 'trans'], required=True, help='Analysis mode')
    
    # Core Options
    parser.add_argument('--window', default='1000000', help='Cis-window size (default 1Mb)')
    parser.add_argument('--permutations', default='1000', help='Number of permutations (for mode=permute)')
    parser.add_argument('--nominal-threshold', default='1.0', help='Significance threshold for nominal pass')
    parser.add_argument('--chunk-index', help='Chunk index (for parallelization)')
    parser.add_argument('--chunk-total', help='Total chunks (for parallelization)')
    
    # Pass-through arguments
    parser.add_argument('extra_args', nargs=argparse.REMAINDER, help='Additional arguments passed to QTLtools (e.g., --seed, --exclude-samples)')
    
    return parser.parse_args()

def main():
    args = get_args()
    check_qtltools()
    
    # Ensure output directory exists
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Construct Command
    cmd = ['QTLtools']
    
    if args.mode == 'trans':
        cmd.append('trans')
    else:
        cmd.append('cis')
        
    cmd.extend(['--vcf', args.vcf])
    cmd.extend(['--bed', args.bed])
    
    if args.cov:
        cmd.extend(['--cov', args.cov])
        
    cmd.extend(['--out', args.out])
    
    # Mode-specific args
    if args.mode == 'permute':
        cmd.extend(['--permute', str(args.permutations)])
        cmd.extend(['--window', str(args.window)])
        
    elif args.mode == 'nominal':
        cmd.extend(['--nominal', str(args.nominal_threshold)])
        cmd.extend(['--window', str(args.window)])
        
    elif args.mode == 'trans':
        cmd.extend(['--permute', str(args.permutations)])
    
    # Parallelization
    if args.chunk_index and args.chunk_total:
        cmd.extend(['--chunk', str(args.chunk_index), str(args.chunk_total)])
        
    # Append any extra arguments
    if args.extra_args:
        print(f"Passing extra arguments to QTLtools: {args.extra_args}")
        cmd.extend(args.extra_args)
        
    # Execution
    print("Running QTLtools...")
    print("Command:", " ".join(cmd))
    
    try:
        subprocess.check_call(cmd)
        print("Analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing QTLtools: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
