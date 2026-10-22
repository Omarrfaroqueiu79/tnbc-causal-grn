"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 1: Sparse ingestion, quality control filtering, and cell cohort stratification for single-cell CRISPRi Perturb-seq data (GSE224986).

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 3]
import os
import gc
import time
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from tqdm.auto import tqdm

sc.settings.verbosity = 3
sc.set_figure_params(dpi=100, color_map='viridis')

# Dynamic root directory resolution
BASE_DIR = r'O:\Thesis\TNBC_Work\TNBC_Final_Dataset'
if not os.path.exists(BASE_DIR):
    BASE_DIR = os.getcwd()

SOURCE_DIR = os.path.join(BASE_DIR, 'source')
TRANS_FILE = os.path.join(SOURCE_DIR, 'GSE224986_MDAMB231_SM_trans_df.h5')
SGRNA_FILE = os.path.join(SOURCE_DIR, 'GSE224986_MDAMB231_SM_sgRNA_df.h5')
OUTPUT_H5AD = os.path.join(BASE_DIR, 'MDA_MB_231_SM_analysis_ready.h5ad')

print(f'Project root: {BASE_DIR}')
print(f'Transcriptome source: {TRANS_FILE} (Exists: {os.path.exists(TRANS_FILE)})')
print(f'sgRNA source: {SGRNA_FILE} (Exists: {os.path.exists(SGRNA_FILE)})')

# %% [Code Cell 5]
start_time = time.time()
print('Reading sgRNA matrix in memory-safe chunks...')

chunk_size_sgrna = 2000
sg_chunks = []

with h5py.File(SGRNA_FILE, 'r') as f:
    sg_ds = f['df']['block0_values']  # Shape: (11975 guides, 95678 cells)
    guides = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in f['df']['axis1'][:]]
    cells_sg = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in f['df']['axis0'][:]]
    
    for i in tqdm(range(0, sg_ds.shape[0], chunk_size_sgrna), desc='Loading sgRNA chunks'):
        chunk = sg_ds[i : i + chunk_size_sgrna, :]
        chunk = np.nan_to_num(chunk, nan=0.0)
        sg_chunks.append(sp.csr_matrix(chunk))

# Stack chunks (guides x cells) and transpose to standard (cells x guides)
sg_sparse_T = sp.vstack(sg_chunks)
del sg_chunks
gc.collect()

sg_sparse = sg_sparse_T.T.tocsr()
del sg_sparse_T
gc.collect()

sg_ram = sg_sparse.data.nbytes / (1024**2)
print(f'sgRNA sparse matrix constructed: {sg_sparse.shape} (RAM: {sg_ram:.2f} MB)')

print('Constructing cell-level metadata table...')
guide_counts_per_cell = np.diff(sg_sparse.indptr)
total_guide_umis = np.array(sg_sparse.sum(axis=1)).flatten()

primary_guides = []
primary_counts = []
guide_categories = []

for i in range(sg_sparse.shape[0]):
    s_idx = sg_sparse.indptr[i]
    e_idx = sg_sparse.indptr[i + 1]
    n_g = e_idx - s_idx
    if n_g == 0:
        primary_guides.append('unassigned')
        primary_counts.append(0.0)
        guide_categories.append('unassigned')
    else:
        data = sg_sparse.data[s_idx:e_idx]
        indices = sg_sparse.indices[s_idx:e_idx]
        top_i = np.argmax(data)
        primary_guides.append(guides[indices[top_i]])
        primary_counts.append(data[top_i])
        guide_categories.append('single_guide' if n_g == 1 else 'multi_guide')

obs_df = pd.DataFrame({
    'cell_id': cells_sg,
    'cell_line': 'MDA-MB-231',
    'screen': 'SM',
    'lane': [c.split('-')[-1] for c in cells_sg],
    'num_guides': guide_counts_per_cell,
    'total_guide_umis': total_guide_umis,
    'primary_guide': primary_guides,
    'primary_guide_count': primary_counts,
    'guide_assignment': guide_categories
}, index=cells_sg)

print('\nGuide Assignment Breakdown:')
display(obs_df['guide_assignment'].value_counts())
obs_ram = obs_df.memory_usage().sum() / (1024**2)
print(f'Metadata table memory: {obs_ram:.2f} MB')

# %% [Code Cell 7]
print('Loading single-cell transcriptome matrix in chunks...')
chunk_size_trans = 5000
tr_chunks = []

with h5py.File(TRANS_FILE, 'r') as f:
    tr_ds = f['df']['block0_values']  # Shape: (58391 genes, 95678 cells)
    genes = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in f['df']['axis1'][:]]
    cells_tr = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in f['df']['axis0'][:]]
    
    for i in tqdm(range(0, tr_ds.shape[0], chunk_size_trans), desc='Loading transcriptome chunks'):
        chunk = tr_ds[i : i + chunk_size_trans, :]
        tr_chunks.append(sp.csr_matrix(chunk))

tr_sparse_T = sp.vstack(tr_chunks)
del tr_chunks
gc.collect()

tr_sparse = tr_sparse_T.T.tocsr()  # (cells x genes)
del tr_sparse_T
gc.collect()

tr_ram = tr_sparse.data.nbytes / (1024**2)
print(f'Transcriptome sparse matrix ready: {tr_sparse.shape} (RAM: {tr_ram:.2f} MB)')

# %% [Code Cell 9]
print('Assembling unified AnnData object...')
adata = ad.AnnData(
    X=tr_sparse,
    obs=obs_df,
)
adata.obs_names = cells_sg
adata.var_names = genes
adata.var_names_make_unique()

adata.layers['counts'] = adata.X.copy()
adata.obsm['sgRNA_counts'] = sg_sparse
adata.uns['sgRNA_names'] = np.array(guides)

# Free intermediate variables and collect garbage
del tr_sparse, sg_sparse, obs_df
gc.collect()

print('AnnData object successfully initialized:')
print(adata)

# %% [Code Cell 11]
print('Calculating QC metrics...')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

print(f'Pre-filter dimensions: {adata.shape}')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f'Post-filter dimensions: {adata.shape}')

gc.collect()

# %% [Code Cell 13]
print(f'Saving finalized AnnData to: {OUTPUT_H5AD}...')
adata.write_h5ad(OUTPUT_H5AD)

file_size_mb = os.path.getsize(OUTPUT_H5AD) / (1024**2)
total_sec = time.time() - start_time
print(f'\nSUCCESS! Analysis-ready dataset finalized in {total_sec:.1f} seconds.')
print(f'Output file: {OUTPUT_H5AD}')
print(f'File size: {file_size_mb:.2f} MB')

# Verification summary
print('\n--- Dataset Summary ---')
print(f'Cells: {adata.n_obs:,}')
print(f'Genes: {adata.n_vars:,}')
print(f'Perturbations (guides): {adata.obsm["sgRNA_counts"].shape[1]:,}')
print(f'Lanes: {sorted(adata.obs["lane"].unique())}')
print(f'Counts layer preserved: {"counts" in adata.layers}')

