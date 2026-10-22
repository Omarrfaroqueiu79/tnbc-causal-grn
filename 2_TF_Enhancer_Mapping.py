"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 2: Multi-omics transcription factor prior mapping against ENCODE Factorbook ChIP-seq peaks across 819 somatic mutation enhancer loci.

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 3]
import os
import gzip
import time
import requests
import numpy as np
import pandas as pd
import anndata as ad
from tqdm.auto import tqdm

BASE_DIR = r'O:\Thesis\TNBC_Work\TNBC_Final_Dataset'
if not os.path.exists(BASE_DIR):
    BASE_DIR = os.getcwd()

H5AD_FILE = os.path.join(BASE_DIR, 'MDA_MB_231_SM_analysis_ready.h5ad')
BED_FILE = os.path.join(BASE_DIR, 'source', 'GSE223812_Somatic_mutation_screen_sgRNAs.bed.gz')
OUTPUT_PRIOR_CSV = os.path.join(BASE_DIR, 'TF_Enhancer_Priors.csv')

print(f'Working directory: {BASE_DIR}')
print(f'Finalized dataset: {H5AD_FILE} (Exists: {os.path.exists(H5AD_FILE)})')
print(f'Enhancer BED file: {BED_FILE} (Exists: {os.path.exists(BED_FILE)})')

# %% [Code Cell 5]
print('Loading finalized AnnData...')
adata = ad.read_h5ad(H5AD_FILE)
guides_in_screen = set(adata.uns['sgRNA_names'])
print(f'Loaded dataset with {adata.n_obs:,} cells and {len(guides_in_screen):,} unique sgRNA guides.')

# Parse guide-to-enhancer coordinates from BED
print('\nParsing enhancer coordinates from somatic mutation screen BED...')
records = []
with gzip.open(BED_FILE, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name_parts = parts[3].split(';')
            guide_id = name_parts[0]
            seq = name_parts[-1]
            enh_part = [p for p in name_parts if '::' in p]
            enh_region = enh_part[0].split('::')[0].replace('>', '') if enh_part else f'{chrom}:{start}-{end}'
            records.append({
                'chrom': chrom, 'start': start, 'end': end,
                'guide_id': guide_id, 'sequence': seq,
                'enhancer_region': enh_region
            })

df_guides = pd.DataFrame(records)
df_matched = df_guides[df_guides['sequence'].isin(guides_in_screen)].copy()
unique_enhancers = df_matched['enhancer_region'].unique()

print(f'Matched {len(df_matched):,} guides to {len(unique_enhancers):,} unique enhancer loci.')
display(df_matched.head(4))

# %% [Code Cell 7]
def get_tfs_for_enhancer(enhancer_str, session=None):
    '''
    Query UCSC ENCODE TF Binding Site clusters (encRegTfbsClustered) for a genomic interval.
    Input format: chr10:102193114-102193961
    '''
    try:
        chrom, coords = enhancer_str.split(':')
        start, end = coords.split('-')
        url = f'https://api.genome.ucsc.edu/getData/track?genome=hg38&track=encRegTfbsClustered&chrom={chrom}&start={start}&end={end}'
        
        s = session or requests
        r = s.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            items = data.get('encRegTfbsClustered', [])
            tf_results = []
            for it in items:
                tf_results.append({
                    'tf_name': it.get('name'),
                    'score': it.get('score'),
                    'binding_start': it.get('chromStart'),
                    'binding_end': it.get('chromEnd')
                })
            return tf_results
        return []
    except Exception:
        return []

# Quick test query
test_enh = unique_enhancers[0]
print(f'Test query for {test_enh}:')
test_tfs = get_tfs_for_enhancer(test_enh)
print(f'Found {len(test_tfs)} TF binding sites: {[t["tf_name"] for t in test_tfs]}')

# %% [Code Cell 9]
print(f'Mapping {len(unique_enhancers):,} enhancers to ENCODE TF clusters...')
prior_edges = []
session = requests.Session()

# Batch query with rate-limiting and progress tracking
for enh in tqdm(unique_enhancers, desc='Querying ENCODE TFBS'):
    chrom, coords = enh.split(':')
    start, end = coords.split('-')
    tfs = get_tfs_for_enhancer(enh, session=session)
    
    for item in tfs:
        prior_edges.append({
            'tf_name': item['tf_name'],
            'enhancer_region': enh,
            'chrom': chrom,
            'enhancer_start': int(start),
            'enhancer_end': int(end),
            'binding_start': item['binding_start'],
            'binding_end': item['binding_end'],
            'confidence_score': item['score'],
            'database': 'ENCODE_ChIP_seq'
        })
    time.sleep(0.05)  # Respect API rate limits

df_priors = pd.DataFrame(prior_edges)
print(f'\nTotal TF -> Enhancer prior edges identified: {len(df_priors):,}')
display(df_priors.head(5))

# %% [Code Cell 11]
df_priors.to_csv(OUTPUT_PRIOR_CSV, index=False)
print(f'Saved TF -> Enhancer prior graph to: {OUTPUT_PRIOR_CSV}')

print('\n--- Summary of Identified Epigenetic Regulators ---')
print(f'Unique Transcription Factors: {df_priors["tf_name"].nunique():,}')
print(f'Unique Enhancers Bound: {df_priors["enhancer_region"].nunique():,} / {len(unique_enhancers):,}')

print('\nTop 15 Most Prevalent Transcription Factors in Perturbed Enhancers:')
display(df_priors['tf_name'].value_counts().head(15))

