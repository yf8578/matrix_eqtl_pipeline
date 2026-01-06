import pandas as pd
import argparse
import sys
import os
import subprocess
import requests
import io

def check_dependencies():
    """Check if bgzip and tabix are available (only needed for QTLtools format)."""
    for cmd in ['bgzip', 'tabix']:
        if subprocess.call(['which', cmd], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            print(f"Warning: {cmd} is not found in PATH. BED file will not be compressed/indexed.")
            return False
    return True

def get_args():
    parser = argparse.ArgumentParser(description="Prepare Expression/Position Data for MatrixEQTL or QTLtools")
    
    # Inputs
    parser.add_argument('--expression', required=True, help='Expression Matrix file (ID in first col)')
    
    # Position Source (Mutually Exclusive-ish)
    parser.add_argument('--positions', help='Local Position file (geneid, chr, left, right)')
    parser.add_argument('--fetch-pos', action='store_true', help='Auto-download positions from BioMart')
    parser.add_argument('--gtf', help='Extract positions from local GTF file')
    parser.add_argument('--dummy-pos', action='store_true', help='Generate dummy positions (1MB spacing)')
    parser.add_argument('--dummy-chr', default='1', help='Chromosome for dummy positions')
    
    # Settings
    parser.add_argument('--id-type', default='ensembl_gene_id', 
                        help='ID type for fetching/GTF: ensembl_gene_id, external_gene_name, entrezgene_id, uniprot_gn_id')
    
    # Output Control
    parser.add_argument('--out', required=True, help='Output filename prefix or full path')
    parser.add_argument('--format', choices=['qtltools', 'matrixeqtl'], default='qtltools',
                        help='Output format. "qtltools": BED.gz with expression. "matrixeqtl": TSV with positions only.')
    
    # Legacy
    parser.add_argument('--strand-ref', help='Optional strand ref file (legacy)')
    
    return parser.parse_args()

def parse_gtf(gtf_file, target_ids, id_type):
    """Parse GTF for coordinates."""
    print(f"Parsing GTF: {gtf_file}...")
    attr_map = {
        'ensembl_gene_id': 'gene_id',
        'external_gene_name': 'gene_name',
        'gene_id': 'gene_id',
        'gene_name': 'gene_name',
        'entrezgene_id': 'db_xref',
        'uniprot_gn_id': 'NA'
    }
    target_attr = attr_map.get(id_type, id_type)
    
    if target_attr == 'NA':
        print(f"Warning: ID type {id_type} not typically supported in GTF parsing.")
        
    results = {}
    target_set = set(target_ids)
    
    opener = open
    if gtf_file.endswith('.gz'):
        import gzip
        opener = gzip.open
        
    try:
        with opener(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.split('\t')
                if len(parts) < 9: continue
                if parts[2] != 'gene': continue
                
                attrs = parts[8]
                # Fast heuristic search
                if target_attr not in attrs: continue
                
                # Extract value: key "value";
                # Find key
                k_idx = attrs.find(target_attr)
                v_start = attrs.find('"', k_idx) + 1
                v_end = attrs.find('"', v_start)
                if v_start == 0 or v_end == -1: continue
                
                val = attrs[v_start:v_end]
                
                if val in target_set:
                    results[val] = [val, parts[0], int(parts[3]), int(parts[4]), parts[6]]
                    
    except Exception as e:
        print(f"Error reading GTF: {e}")
        return pd.DataFrame()
        
    if not results: return pd.DataFrame()
    return pd.DataFrame.from_dict(results, orient='index', columns=['geneid', 'chr', 'left', 'right', 'strand'])

def fetch_biomart_positions(ids, id_type='ensembl_gene_id'):
    print(f"Querying BioMart for {len(ids)} IDs ({id_type})...")
    
    # Map to BioMart Attributes
    # List: http://www.ensembl.org/biomart/martview
    mart_attr_map = {
        'ensembl_gene_id': 'ensembl_gene_id',
        'ensembl': 'ensembl_gene_id',
        'external_gene_name': 'external_gene_name',
        'symbol': 'external_gene_name',
        'entrezgene_id': 'entrezgene_id',
        'entrez': 'entrezgene_id',
        'uniprot_gn_id': 'uniprot_gn_id', # UniProt Gene Name
        'uniprot': 'uniprot_gn_id'
        # Note: UniProt Accession is 'uniprot_gn_symbol' or 'uniprotswissprot' etc.
    }
    query_attr = mart_attr_map.get(id_type, id_type)
    
    chunk_size = 200
    unique_ids = list(set(ids))
    all_res = []
    
    xml_template = """<!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "{filter_name}" value = "{val}"/>
            <Attribute name = "{filter_name}" />
            <Attribute name = "chromosome_name" />
            <Attribute name = "start_position" />
            <Attribute name = "end_position" />
            <Attribute name = "strand" />
        </Dataset>
    </Query>"""
    
    for i in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[i:i+chunk_size]
        vals = ",".join([str(x).strip() for x in chunk])
        query = xml_template.format(filter_name=query_attr, val=vals)
        
        try:
            r = requests.post("http://www.ensembl.org/biomart/martservice", data={'query': query})
            r.raise_for_status()
            if not r.text.strip(): continue
            
            # Read
            df = pd.read_csv(io.StringIO(r.text), sep='\t', header=None, dtype=str)
            # Expect 5 cols
            if df.shape[1] >= 5:
                df = df.iloc[:, :5]
                df.columns = ['geneid', 'chr', 'left', 'right', 'strand']
                all_res.append(df)
        except Exception as e:
            print(f"Warning: Batch failed: {e}")
            
    if not all_res: return pd.DataFrame()
    final = pd.concat(all_res)
    
    # Filter standard chromosomes
    valid = set([str(x) for x in range(1,23)] + ['X','Y','MT','M'])
    final = final[final['chr'].isin(valid)]
    
    # Dedup
    final = final.drop_duplicates(subset=['geneid'])
    
    return final

def main():
    args = get_args()
    has_tools = check_dependencies()
    
    # 1. Read Expression
    print(f"Reading Expression: {args.expression}")
    expr = pd.read_csv(args.expression, sep='\t')
    # Rename ID col
    expr.rename(columns={expr.columns[0]: 'geneid'}, inplace=True)
    expr['geneid'] = expr['geneid'].astype(str)
    
    ids = expr['geneid'].tolist()
    
    # 2. Get Positions
    pos_df = pd.DataFrame()
    
    if args.dummy_pos:
        print("Generating dummy positions...")
        starts = [1000 + i*1000 for i in range(len(ids))]
        pos_df = pd.DataFrame({
            'geneid': ids,
            'chr': [args.dummy_chr]*len(ids),
            'left': starts,
            'right': [s+100 for s in starts],
            'strand': ['.']*len(ids)
        })
    elif args.fetch_pos:
        pos_df = fetch_biomart_positions(ids, args.id_type)
    elif args.gtf:
        if not os.path.exists(args.gtf):
             print(f"Error: GTF not found: {args.gtf}")
             sys.exit(1)
        pos_df = parse_gtf(args.gtf, ids, args.id_type)
    elif args.positions:
        print(f"Reading positions file: {args.positions}")
        pos_df = pd.read_csv(args.positions, sep='\t')
        pos_df.columns = [c.lower() for c in pos_df.columns]
        # Map cols
        if 'features' in pos_df.columns: pos_df.rename(columns={'features': 'geneid'}, inplace=True)
        if 'start' in pos_df.columns: pos_df.rename(columns={'start': 'left'}, inplace=True)
        if 'end' in pos_df.columns: pos_df.rename(columns={'end': 'right'}, inplace=True)
        if 's1' in pos_df.columns: pos_df.rename(columns={'s1': 'left'}, inplace=True)
        if 's2' in pos_df.columns: pos_df.rename(columns={'s2': 'right'}, inplace=True)
    else:
        print("Error: No position source specified.")
        sys.exit(1)
        
    if pos_df.empty:
        print("Error: No positions found/matched.")
        sys.exit(1)
        
    # Standardize Pos DF
    pos_df['geneid'] = pos_df['geneid'].astype(str)
    if 'strand' not in pos_df.columns: pos_df['strand'] = '.'
    
    # 3. Merge
    print("Merging features...")
    # Inner join to keep only matching genes
    merged = pd.merge(expr, pos_df, on='geneid', how='inner')
    print(f"  Matched {len(merged)} features.")
    
    # 4. Output
    if args.format == 'matrixeqtl':
        out_file = args.out
        if not out_file.endswith('.txt') and not out_file.endswith('.tsv'):
            out_file += '.txt'
            
        print(f"Writing MatrixEQTL Position file to: {out_file}")
        # Columns: geneid, chr, left, right
        final = merged[['geneid', 'chr', 'left', 'right']]
        final.to_csv(out_file, sep='\t', index=False)
        print("Done.")
        
    else: # qtltools
        print("Preparing BED format...")
        # Columns: #Chr, start, end, pid, gid, strand, samples...
        bed = pd.DataFrame()
        bed['#Chr'] = merged['chr'].astype(str)
        # BED is 0-based start
        bed['start'] = merged['left'].astype(int) - 1
        bed['end'] = merged['right'].astype(int)
        bed['pid'] = merged['geneid']
        bed['gid'] = merged['geneid']
        
        # Strand Format
        def fix_strand(x):
            if x in ['1', '+']: return '+'
            if x in ['-1', '-']: return '-'
            return '.'
        bed['strand'] = merged['strand'].apply(fix_strand)
        
        # Samples
        samples = [c for c in expr.columns if c != 'geneid']
        bed = pd.concat([bed, merged[samples]], axis=1)
        
        # Sort
        bed.sort_values(by=['#Chr', 'start'], ascending=[True, True], inplace=True)
        
        out_bed = args.out
        if not out_bed.endswith('.bed'): out_bed += '.bed'
        
        print(f"Writing BED: {out_bed}")
        bed.to_csv(out_bed, sep='\t', index=False, float_format='%.5g')
        
        if has_tools:
            print("Compressing & Indexing...")
            subprocess.check_call(f"bgzip -f {out_bed}", shell=True)
            subprocess.check_call(f"tabix -p bed {out_bed}.gz", shell=True)
            print(f"Output: {out_bed}.gz")
        else:
            print(f"Output: {out_bed} (Uncompressed)")

if __name__ == "__main__":
    main()
