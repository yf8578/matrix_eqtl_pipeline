import argparse
import requests
import pandas as pd
import io
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Fetch Gene Positions from BioMart (Run Locally)")
    parser.add_argument('--ids', help='File containing list of IDs (one per line) OR comma-separated string')
    parser.add_argument('--id-type', default='ensembl_gene_id', help='ID Type: ensembl_gene_id, external_gene_name, etc.')
    parser.add_argument('--genome', default='hg38', choices=['hg38', 'hg19'], help='Genome Build (hg38=Ensembl, hg19=Grch37)')
    parser.add_argument('--out', default='gene_positions.tsv', help='Output filename')
    parser.add_argument('--all-genes', action='store_true', help='Download ALL genes (no ID filtering)')
    return parser.parse_args()

def fetch_biomart(ids=None, id_type='ensembl_gene_id', genome='hg38'):
    print(f"Connecting to BioMart ({genome})...")
    
    # URL Selection
    if genome == 'hg19':
        url = "http://grch37.ensembl.org/biomart/martservice"
        dataset = "hsapiens_gene_ensembl"
    else:
        url = "http://www.ensembl.org/biomart/martservice"
        dataset = "hsapiens_gene_ensembl"

    # Attributes
    id_attr = id_type
    
    xml = f"""
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "{dataset}" interface = "default" >
            {f'<Filter name = "{id_attr}" value = "{",".join(ids)}"/>' if ids else ''}
            <Attribute name = "{id_attr}" />
            <Attribute name = "chromosome_name" />
            <Attribute name = "start_position" />
            <Attribute name = "end_position" />
            <Attribute name = "strand" />
        </Dataset>
    </Query>
    """
    
    try:
        r = requests.post(url, data={'query': xml}, stream=True)
        r.raise_for_status()
        
        # Determine columns based on header
        # Usually: ID, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand
        df = pd.read_csv(io.StringIO(r.text), sep='\t', dtype=str)
        
        # Standardize Columns
        # Rename to: geneid, chr, left, right, strand
        # The column names from BioMart are verbose, e.g. "Gene stable ID", "Chromosome/scaffold name"
        # We rename by position to be safe (assuming attr order matches query)
        # XML requests usually return in order requested
        
        new_cols = ['geneid', 'chr', 'left', 'right', 'strand']
        if len(df.columns) == 5:
            df.columns = new_cols
        else:
            # Fallback mapping if BioMart adds extra columns
            print("Warning: unexpected column count. Using raw columns.")
            print(df.columns)
            
        return df
        
    except Exception as e:
        print(f"Error: {e}")
        return pd.DataFrame()

def main():
    args = get_args()
    
    id_list = []
    if args.ids:
        if ',' in args.ids:
            id_list = [x.strip() for x in args.ids.split(',')]
        elif os.path.exists(args.ids):
            with open(args.ids) as f:
                id_list = [line.strip().split()[0] for line in f if line.strip()]
        else:
            print("Error: --ids should be a comma-separated string or a valid file path.")
            sys.exit(1)
            
    if not id_list and not args.all_genes:
        print("Error: Must specify --ids or --all-genes")
        sys.exit(1)
        
    print(f"Fetching positions for {len(id_list) if id_list else 'ALL'} genes...")
    
    # Fetch
    # If list is too long, we might need chunking, but for simple user utility let's try one go or chunk simple
    # The XML post method handles reasonably large bodies.
    
    if args.all_genes:
         df = fetch_biomart(None, args.id_type, args.genome)
    else:
         # Chunk 500
         chunk_size = 500
         dfs = []
         unique_ids = list(set(id_list))
         for i in range(0, len(unique_ids), chunk_size):
             print(f"  Batch {i}-{i+chunk_size}...")
             sub_ids = unique_ids[i:i+chunk_size]
             dfs.append(fetch_biomart(sub_ids, args.id_type, args.genome))
         df = pd.concat(dfs)

    if df.empty:
        print("No data retrieved.")
        sys.exit(1)
        
    # Clean
    valid_chrs = set([str(c) for c in range(1, 23)] + ['X', 'Y', 'MT', 'M'])
    df = df[df['chr'].isin(valid_chrs)]
    
    # Save
    print(f"Saving {len(df)} records to {args.out}...")
    df.to_csv(args.out, sep='\t', index=False)
    print("Done! Upload this file to your cluster and use with --positions.")

if __name__ == "__main__":
    import os
    main()
