import pandas as pd
import argparse
import sys
import os
import subprocess

def check_dependencies():
    """Check if bgzip and tabix are available."""
    for cmd in ['bgzip', 'tabix']:
        if subprocess.call(['which', cmd], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"Warning: {cmd} is not found. BED files will not be compressed/indexed.")

def get_args():
    parser = argparse.ArgumentParser(description="Map Features to Genomic Locations and Prepare Inputs")
    
    # Inputs
    parser.add_argument('--matrix', required=True, help='Feature Matrix (Rows=Features, Cols=Samples). First col must be ID.')
    parser.add_argument('--master-loc', required=True, help='Master Gene Location File (geneid, chr, left, right, [strand])')
    
    # ID Mapping Options
    parser.add_argument('--input-id', default='ensembl', choices=['ensembl', 'symbol', 'entrez', 'uniprot'], help='ID type in your Matrix')
    parser.add_argument('--ref-id', default='ensembl', choices=['ensembl', 'symbol', 'entrez'], help='ID type in Master Location File')
    parser.add_argument('--auto-map', action='store_true', help='Automatically fetch mapping from BioMart if IDs differ')
    
    # Output Control
    parser.add_argument('--format', choices=['matrixeqtl', 'qtltools'], required=True, help='Output format')
    parser.add_argument('--out', required=True, help='Output file path')
    
    return parser.parse_args()

def get_biomart_attr_name(friendly_name):
    """Map friendly names to Ensembl BioMart attribute names."""
    mapping = {
        'ensembl': 'ensembl_gene_id',
        'symbol': 'external_gene_name',
        'entrez': 'entrezgene_id',
        'uniprot': 'uniprot_gn_id' # Swissprot/TrEMBL ID
    }
    return mapping.get(friendly_name, friendly_name)

def main():
    args = get_args()
    check_dependencies()
    
    # 1. Read Feature Matrix
    print(f"Reading Feature Matrix: {args.matrix}")
    expr = pd.read_csv(args.matrix, sep='\t')
    feat_col = expr.columns[0]
    expr.rename(columns={feat_col: 'feature_id'}, inplace=True)
    # Ensure ID is string
    expr['feature_id'] = expr['feature_id'].astype(str)
    
    features = expr[['feature_id']].drop_duplicates()
    print(f"  Found {len(features)} unique features.")

    # 2. Read Master Location File
    print(f"Reading Master Locations: {args.master_loc}")
    master = pd.read_csv(args.master_loc, sep='\t')
    master.columns = [c.lower() for c in master.columns]
    
    # Standardize Master Columns
    col_map = {'features': 'ref_id', 'geneid': 'ref_id', 'gene_id':'ref_id', 'start':'left', 'end':'right', 'chromosome':'chr'}
    master.rename(columns=col_map, inplace=True)
    master['ref_id'] = master['ref_id'].astype(str)
    
    if 'ref_id' not in master.columns:
        print("Error: Master location file must have a 'geneid' (or features/gene_id) column.")
        sys.exit(1)

    # 3. Perform Mapping
    mapped_data = None
    
    # If IDs are same type, merge directly
    if args.input_id == args.ref_id:
        print(f"ID types match ({args.input_id}). Merging directly...")
        # feature_id IS ref_id
        features['ref_id'] = features['feature_id']
        mapped_data = pd.merge(features, master, on='ref_id', how='inner')
        
    else:
        # IDs match, need mapping
        if not args.auto_map:
            print("Error: ID types differ but --auto-map is not enabled.")
            print(f"  Matrix: {args.input_id}, Master: {args.ref_id}")
            sys.exit(1)
            
        print(f"Need mapping: {args.input_id} -> {args.ref_id}")
        
        # Prepare Auto-Fetch
        source_attr = get_biomart_attr_name(args.input_id)
        target_attr = get_biomart_attr_name(args.ref_id)
        
        # Check if we also need to fetch strand (if missing in Master)
        extra_attrs = []
        fetch_strand = False
        if 'strand' not in master.columns:
            print("  Master file lacks 'strand'. Will attempt to fetch 'strand' from BioMart.")
            extra_attrs.append('strand')
            fetch_strand = True
        
        map_file = args.out + ".mapping_temp.txt"
        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script = os.path.join(script_dir, "fetch_biomart_mapping.R")
        
        print(f">>> Calling BioMart to fetch mapping table (plus {extra_attrs})...")
        cmd = ['Rscript', r_script, map_file, source_attr, target_attr]
        # dataset name arg (optional in R script, but positional at index 4)
        # We need to pass dataset name if we pass extras
        cmd.append('hsapiens_gene_ensembl') # Default dataset
        cmd.extend(extra_attrs)
        
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            print("Error: Failed to fetch mapping from BioMart.")
            if os.path.exists(map_file): os.remove(map_file)
            sys.exit(1)
            
        # Load Mapping
        print("Loading mapping table...")
        id_map = pd.read_csv(map_file, sep='\t')
        
        # Cleanup
        if os.path.exists(map_file): os.remove(map_file)

        # Standardize Map Columns
        # Typically: input_id_col, ref_id_col, [strand]
        # We know source_attr -> feature_id, target_attr -> ref_id
        # Note: BioMart column names match attributes
        rename_map = {source_attr: 'feature_id', target_attr: 'ref_id'}
        if fetch_strand: rename_map['strand'] = 'strand_fetched'
        
        id_map.rename(columns=rename_map, inplace=True)
        id_map['feature_id'] = id_map['feature_id'].astype(str)
        id_map['ref_id'] = id_map['ref_id'].astype(str)
        
        # Merge: Feature(Matrix) -> Map -> MasterLoc
        step1 = pd.merge(features, id_map, on='feature_id', how='inner')
        print(f"  Mapped {len(step1)} features to reference IDs.")
        
        mapped_data = pd.merge(step1, master, on='ref_id', how='inner')
        
        # If we fetched strand, use it
        if fetch_strand and 'strand_fetched' in mapped_data.columns:
            mapped_data['strand'] = mapped_data['strand_fetched']

    print(f"  Final matched locations: {len(mapped_data)}")
    
    if len(mapped_data) == 0:
        print("Error: No features matched coordinates. Check your IDs.")
        sys.exit(1)

    # 4. Generate Output (mapped_data has feature_id, ref_id, chr, left, right...)
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir): os.makedirs(out_dir)

    if args.format == 'matrixeqtl':
        # Output: ID, Chr, Left, Right
        # ID is the ORIGINAL feature_id (e.g. Protein ID)
        out_df = mapped_data[['feature_id', 'chr', 'left', 'right']].copy()
        out_df.columns = ['id', 'chr', 's1', 's2']
        # Remove duplicates if one protein maps to multiple genes (rare but possible? or vice versa)
        # Usually we want unique feature_ids. If a protein maps to 2 genes, MatrixEQTL might fail if IDs duplicate?
        # MatrixEQTL location file ID must be unique.
        out_df.drop_duplicates(subset=['id'], inplace=True)
        
        print(f"Writing MatrixEQTL location file: {args.out}")
        out_df.to_csv(args.out, sep='\t', index=False)
        
    elif args.format == 'qtltools':
        check_dependencies()
        
        # Merge Coordinates back to Expression Data
        # mapped_data: [feature_id, ref_id, chr, left, right, strand...]
        full_data = pd.merge(mapped_data, expr, on='feature_id', how='inner')
        
        # Deduplicate features (if mapping yielded multi-hits, pick one?)
        # For safety, let's just drop duplicates on feature_id
        full_data.drop_duplicates(subset=['feature_id'], inplace=True)
        
        # Build BED
        bed = pd.DataFrame()
        bed['#Chr'] = full_data['chr'].astype(str)
        bed['start'] = full_data['left'].astype(int)
        bed['end'] = full_data['right'].astype(int)
        bed['pid'] = full_data['feature_id'] # Phenotype ID = Feature ID (e.g. UniProt)
        bed['gid'] = full_data['ref_id']     # Group ID = Gene ID (e.g. Ensembl)
        
        # Strand
        if 'strand' in full_data.columns:
            def clean_strand(x):
                x = str(x).strip()
                if x in ['1','+']: return '+'
                if x in ['-1','-']: return '-'
                return '.'
            bed['strand'] = full_data['strand'].apply(clean_strand)
        else:
            bed['strand'] = '.'
            
        # Samples
        sample_cols = [c for c in expr.columns if c != 'feature_id']
        bed = pd.concat([bed, full_data[sample_cols]], axis=1)
        
        # Sort
        bed.sort_values(by=['#Chr', 'start'], inplace=True)
        
        # Write
        raw_out = args.out
        if raw_out.endswith('.gz'): raw_out = raw_out[:-3]
        
        print(f"Writing BED: {raw_out}")
        bed.to_csv(raw_out, sep='\t', index=False, float_format='%.5g')
        
        # Compress
        print("Compressing & Indexing...")
        if os.path.exists(raw_out + ".gz"): os.remove(raw_out + ".gz")
        subprocess.check_call(f"bgzip {raw_out}", shell=True)
        subprocess.check_call(f"tabix -p bed {raw_out}.gz", shell=True)
        print(f"Success: {raw_out}.gz")

if __name__ == "__main__":
    main()
