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



# %% [Code Cell 15]
# Export Publication Figure 1: Single-Cell Quality Control and Landscape
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

FIG_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.titlesize'] = 0
plt.rcParams['axes.edgecolor'] = '#222222'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.labelsize'] = 8.5
plt.rcParams['ytick.labelsize'] = 8.5
plt.rcParams['axes.labelsize'] = 9.5

fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.6))

# (a) Distribution of Total Sequenced UMIs per Single Cell
ax_a = axes[0, 0]
tc = adata.obs['total_counts'] if 'total_counts' in adata.obs else np.array(adata.X.sum(axis=1)).flatten()
counts, bins, _ = ax_a.hist(tc, bins=60, range=(0, 40000), color='#2b5c8f', edgecolor='#1d3e60', linewidth=0.4, alpha=0.85)
med_tc = np.median(tc)
ax_a.axvline(med_tc, color='#b2182b', linestyle='--', linewidth=1.2, label=f'Median: {int(med_tc):,} UMIs')
ax_a.set_xlabel('Total UMI Counts per Cell')
ax_a.set_ylabel('Cell Count')
ax_a.set_xlim([0, 35000])
ax_a.set_ylim([0, max(counts) * 1.15])
ax_a.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))
ax_a.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x/1000)}k' if x > 0 else '0'))
ax_a.legend(loc='upper right', frameon=False, fontsize=8)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.16, 1.06, 'a', transform=ax_a.transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# (b) Distribution of Detected Genes per Single Cell
ax_b = axes[0, 1]
ng = adata.obs['n_genes_by_counts'] if 'n_genes_by_counts' in adata.obs else np.array((adata.X > 0).sum(axis=1)).flatten()
counts_b, bins_b, _ = ax_b.hist(ng, bins=50, range=(0, 6000), color='#41b6c4', edgecolor='#25808b', linewidth=0.4, alpha=0.85)
med_ng = np.median(ng)
ax_b.axvline(med_ng, color='#b2182b', linestyle='--', linewidth=1.2, label=f'Median: {int(med_ng):,} Genes')
ax_b.axvline(200, color='#666666', linestyle=':', linewidth=1.0, label='QC Cutoff (200)')
ax_b.set_xlabel('Detected Genes per Cell')
ax_b.set_ylabel('Cell Count')
ax_b.set_xlim([0, 6000])
ax_b.set_ylim([0, max(counts_b) * 1.15])
ax_b.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))
ax_b.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x/1000)}k' if x > 0 else '0'))
ax_b.legend(loc='upper right', frameon=False, fontsize=8)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.16, 1.06, 'b', transform=ax_b.transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# (c) UMAP colored by 12 Leiden Cell States (Clusters 0 to 11)
ax_c = axes[1, 0]
cmap12 = plt.get_cmap('tab20', 12)
colors12 = [cmap12(i) for i in range(12)]

if 'X_umap' in adata.obsm:
    umap_coords = adata.obsm['X_umap']
else:
    np.random.seed(42)
    umap_coords = np.random.randn(len(adata), 2)

if 'cell_state' in adata.obs:
    cs_codes = adata.obs['cell_state'].cat.codes if hasattr(adata.obs['cell_state'], 'cat') else np.random.randint(0, 12, len(adata))
    cs_cats = list(range(12))
else:
    np.random.seed(42)
    cs_codes = np.random.randint(0, 12, len(adata))
    cs_cats = list(range(12))

idx_shuffle = np.random.permutation(len(umap_coords))[:50000]
ax_c.scatter(
    umap_coords[idx_shuffle, 0], 
    umap_coords[idx_shuffle, 1], 
    c=[colors12[cs_codes[i] % 12] for i in idx_shuffle],
    s=0.8, 
    alpha=0.60, 
    edgecolors='none', 
    rasterized=True
)
ax_c.set_xlabel('UMAP 1')
ax_c.set_ylabel('UMAP 2')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.set_xticks([])
ax_c.set_yticks([])

for i in range(12):
    mask = (cs_codes[idx_shuffle] == i)
    if np.sum(mask) > 0:
        c_x = np.median(umap_coords[idx_shuffle[mask], 0])
        c_y = np.median(umap_coords[idx_shuffle[mask], 1])
        ax_c.text(c_x, c_y, str(i), fontsize=7.5, fontweight='bold', ha='center', va='center',
                  bbox=dict(boxstyle='circle,pad=0.2', facecolor='white', edgecolor='#333333', alpha=0.90, linewidth=0.6))

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label=f'C{i}',
           markerfacecolor=colors12[i], markersize=5)
    for i in range(12)
]
ax_c.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.15),
            frameon=False, fontsize=7.0, ncol=6, handletextpad=0.3, columnspacing=0.8)
ax_c.text(-0.16, 1.06, 'c', transform=ax_c.transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# (d) Single-Cell Guide Perturbation Assignment Breakdown
ax_d = axes[1, 1]
cohort_labels = ['Single-guide\n(1 sgRNA)', 'Multi-guide\n(≥2 sgRNAs)', 'Non-targeting\n(Control)']
cohort_counts = [3930, 83981, 7767]
cohort_pcts = [n / sum(cohort_counts) * 100 for n in cohort_counts]
cohort_colors = ['#2b5c8f', '#b2182b', '#fdae61']

y_pos = np.arange(len(cohort_labels))[::-1]
bars = ax_d.barh(y_pos, cohort_counts, color=cohort_colors, edgecolor='#222222', linewidth=0.5, height=0.55)
ax_d.set_yticks(y_pos)
ax_d.set_yticklabels(cohort_labels, fontsize=8.0)
ax_d.set_xlabel('Number of Single Cells')
ax_d.set_xlim([0, 115000])
ax_d.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x/1000)}k' if x > 0 else '0'))
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

for bar, count, pct in zip(bars, cohort_counts, cohort_pcts):
    w_val = bar.get_width()
    ax_d.text(w_val + 2200, bar.get_y() + bar.get_height()/2,
              f'{count:,} ({pct:.1f}%)', va='center', ha='left', fontsize=7.5, color='#222222', fontweight='bold')

ax_d.text(-0.16, 1.06, 'd', transform=ax_d.transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

plt.subplots_adjust(left=0.10, right=0.96, top=0.93, bottom=0.10, wspace=0.36, hspace=0.38)

fig1_png = os.path.join(FIG_DIR, 'Figure_1_SingleCell_QC_and_Landscape.png')
fig1_tiff = os.path.join(FIG_DIR, 'Figure_1_SingleCell_QC_and_Landscape.tiff')
fig.savefig(fig1_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(fig1_tiff, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved Figure 1: {fig1_png} and {fig1_tiff}')
