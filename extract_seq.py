"""
SCRIPT PARAMETERS
path - csv file with the columns
    'gene_id' - Gene id or other label
    'start' - Start basepair
    'end - End basepair
    'strand' - Direction of the sequence
    'chromosome' - The chromosome the gene is on
"""
import sys
import pandas as pd
from tqdm import tqdm
import requests

path = sys.argv[1]
url = 'https://rest.ensembl.org/sequence/region/drosophila_melanogaster/{}:{}..{}'
df = pd.read_csv(path)
out_dict = {
    'gene_id': [],
    'sequence': []
}

for i, row in tqdm(df.iterrows(), total=df.shape[0]):
    start_bp = row['start']
    end_bp = row['stop']
    seq_url = url.format(row['chromosome'], start_bp, end_bp-1)
    r = requests.get(seq_url, headers={ "Content-Type" : "text/x-fasta"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    r = r.text
    if row['strand'] == '-':
        r = r[::-1]
    start = r.index('\n') + 1
    end = r.index('\n') + 1 + (end_bp-start_bp)
    promoter_seq = r[start:end]
    
    out_dict['gene_id'].append(row['gene_id'])
    out_dict['sequence'].append(promoter_seq)

out_df = pd.DataFrame(out_dict)
out_df.to_csv(path[:-4]+' sequences.csv')