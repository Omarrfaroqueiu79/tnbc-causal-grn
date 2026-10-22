"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 3: Normalization, HVG dispersion binning, Leiden phenotypic clustering, and vectorized two-sample Welch t-testing for perturbation responses.

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 3]
import os
import gzip
import time
import gc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

sc.settings.verbosity = 3
sc.set_figure_params(dpi=100, frameon=False, color_map='viridis')

BASE_DIR = r'O:\Thesis\TNBC_Work\TNBC_Final_Dataset'
if not os.path.exists(BASE_DIR):
    BASE_DIR = os.getcwd()

H5AD_FILE = os.path.join(BASE_DIR, 'MDA_MB_231_SM_analysis_ready.h5ad')
BED_FILE = os.path.join(BASE_DIR, 'source', 'GSE223812_Somatic_mutation_screen_sgRNAs.bed.gz')
OUTPUT_ANNOTATED_H5AD = os.path.join(BASE_DIR, 'MDA_MB_231_SM_annotated.h5ad')
OUTPUT_RESPONSES_CSV = os.path.join(BASE_DIR, 'Enhancer_TargetGene_Responses.csv')

print(f'Working directory: {BASE_DIR}')
print(f'Input dataset: {H5AD_FILE} (Exists: {os.path.exists(H5AD_FILE)})')
print(f'Enhancer BED file: {BED_FILE} (Exists: {os.path.exists(BED_FILE)})')

# %% [Code Cell 5]
print('Loading analysis-ready AnnData...')
adata = ad.read_h5ad(H5AD_FILE)
print(f'Loaded {adata.n_obs:,} cells and {adata.n_vars:,} genes.')

# Build guide-to-enhancer dictionary from BED
print('\nMapping guides to enhancer loci...')
guide_to_enh = {}
with gzip.open(BED_FILE, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            chrom, start, end = parts[0], parts[1], parts[2]
            name_parts = parts[3].split(';')
            seq = name_parts[-1]
            enh_part = [p for p in name_parts if '::' in p]
            enh_region = enh_part[0].split('::')[0].replace('>', '') if enh_part else f'{chrom}:{start}-{end}'
            guide_to_enh[seq] = enh_region

# Annotate each cell with its target enhancer
adata.obs['target_enhancer'] = adata.obs['primary_guide'].map(guide_to_enh).fillna('unassigned')

print(f'Mapped {len(guide_to_enh):,} guides to enhancer loci.')
print('Top perturbed enhancer loci by cell count:')
display(adata.obs[adata.obs['guide_assignment'] == 'single_guide']['target_enhancer'].value_counts().head(5))

# %% [Code Cell 7]
print('1. Normalizing expression to 10,000 counts per cell...')
sc.pp.normalize_total(adata, target_sum=1e4)

print('2. Log1p transformation...')
sc.pp.log1p(adata)

print('3. Selecting Highly Variable Genes (HVGs)...')
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False)
print(f'Identified {adata.var["highly_variable"].sum():,} highly variable genes.')

print('4. Running PCA...')
sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')

print('5. Computing neighborhood graph (15 neighbors, 30 PCs)...')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

print('6. Computing UMAP projection...')
sc.tl.umap(adata)

print('7. Performing Leiden community clustering...')
sc.tl.leiden(adata, resolution=0.5, key_added='cell_state', flavor='igraph', n_iterations=2, directed=False)
print(f'Identified {adata.obs["cell_state"].nunique()} distinct TNBC cell states.')

# Visualize UMAP
sc.pl.umap(adata, color=['cell_state', 'guide_assignment', 'lane'], ncols=3, wspace=0.3)

# %% [Code Cell 9]
print('Mapping enhancer perturbation transcriptional responses...')

# Define control/unperturbed cell pool
control_mask = adata.obs['guide_assignment'] == 'unassigned'
single_mask = adata.obs['guide_assignment'] == 'single_guide'

control_cells_idx = np.where(control_mask)[0]
single_cells = adata[single_mask].copy()

# Count cells per targeted enhancer in single-guide cohort
enh_counts = single_cells.obs['target_enhancer'].value_counts()
testable_enhancers = enh_counts[enh_counts >= 3].index.tolist()  # Enhancers with >= 3 cells
testable_enhancers = [e for e in testable_enhancers if e != 'unassigned']

print(f'Total testable enhancers (>=3 single-guide cells): {len(testable_enhancers):,}')

# Pre-compute control mean and standard error for fast vectorized testing
ctrl_X = adata.X[control_cells_idx, :]
ctrl_mean = np.array(ctrl_X.mean(axis=0)).flatten()
ctrl_var = np.array(ctrl_X.power(2).mean(axis=0) - np.square(ctrl_mean)).flatten()
n_ctrl = len(control_cells_idx)

response_records = []
gene_names = np.array(adata.var_names)

# Test each enhancer perturbation for differential response
for enh in tqdm(testable_enhancers, desc='Differential Response Testing'):
    enh_cell_indices = np.where(adata.obs['target_enhancer'] == enh)[0]
    n_pert = len(enh_cell_indices)
    
    pert_X = adata.X[enh_cell_indices, :]
    pert_mean = np.array(pert_X.mean(axis=0)).flatten()
    pert_var = np.array(pert_X.power(2).mean(axis=0) - np.square(pert_mean)).flatten()
    
    # Log2 fold change (difference of log-transformed means)
    log2_fc = pert_mean - ctrl_mean
    
    # Two-sample Welch t-statistic
    denom = np.sqrt(pert_var / n_pert + ctrl_var / n_ctrl + 1e-9)
    t_stat = log2_fc / denom
    
    # Filter candidate responding genes (|log2FC| >= 0.25 and |t_stat| >= 2.5)
    sig_idx = np.where((np.abs(log2_fc) >= 0.25) & (np.abs(t_stat) >= 2.5))[0]
    
    for g_idx in sig_idx:
        response_records.append({
            'enhancer_region': enh,
            'target_gene': gene_names[g_idx],
            'log2_fc': float(log2_fc[g_idx]),
            't_stat': float(t_stat[g_idx]),
            'n_pert_cells': n_pert,
            'effect_direction': 'downregulated' if log2_fc[g_idx] < 0 else 'upregulated'
        })

df_responses = pd.DataFrame(response_records)
print(f'\nIdentified {len(df_responses):,} significant Enhancer -> Target Gene response edges.')
display(df_responses.head(5))

# %% [Code Cell 11]
# 1. Save response edge table
df_responses.to_csv(OUTPUT_RESPONSES_CSV, index=False)
print(f'Saved candidate response edges to: {OUTPUT_RESPONSES_CSV}')

# 2. Save annotated AnnData
print(f'Saving annotated AnnData to: {OUTPUT_ANNOTATED_H5AD}...')
adata.write_h5ad(OUTPUT_ANNOTATED_H5AD)

# 3. Export Publication-Standard Figure: Enhancer Perturbation Responses (Nature/Cell style)
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.titlesize'] = 0
plt.rcParams['axes.edgecolor'] = '#222222'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 9

FIG_DIR = os.path.join(BASE_DIR, 'manuscript_figures')
os.makedirs(FIG_DIR, exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))

# (a) Transcriptional Effect Size (log2FC)
sns.histplot(df_responses['log2_fc'], bins=50, kde=True, ax=axes[0], color='#d6604d', edgecolor='#222222', linewidth=0.5)
axes[0].axvline(0, color='#222222', linestyle='--', linewidth=0.8)
axes[0].set_xlabel('Transcriptional Effect Size (log2 Fold Change)')
axes[0].set_ylabel('Response Edge Count')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].text(-0.16, 1.06, 'a', transform=axes[0].transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# (b) Volcano-style Response Landscape (log2FC vs Welch t-statistic)
sample_resp = df_responses.sample(min(25000, len(df_responses)), random_state=42)
axes[1].scatter(sample_resp['log2_fc'], sample_resp['t_stat'], s=2.5, color='#4575b4', alpha=0.35, rasterized=True)
axes[1].axhline(0, color='#666666', linestyle='--', linewidth=0.7)
axes[1].axvline(0, color='#666666', linestyle='--', linewidth=0.7)
axes[1].set_xlabel('log2 Fold Change')
axes[1].set_ylabel('Welch t-statistic')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].text(-0.14, 1.06, 'b', transform=axes[1].transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

plt.tight_layout()
fig_png = os.path.join(FIG_DIR, 'Figure_3_Enhancer_Perturbation_Responses.png')
fig_tiff = os.path.join(FIG_DIR, 'Figure_3_Enhancer_Perturbation_Responses.tiff')
fig.savefig(fig_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(fig_tiff, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Exported Figure 3: {fig_png} & {fig_tiff} (300 DPI, with panel labels a & b)')

# Summary
print('\n--- Summary of Phase II (Part B) ---')
print(f'Total Analyzed Enhancers: {df_responses["enhancer_region"].nunique():,}')
print(f'Total Responding Target Genes: {df_responses["target_gene"].nunique():,}')
print(f'Mean Response Genes per Enhancer: {len(df_responses) / df_responses["enhancer_region"].nunique():.1f}')
print(f'Saved annotated dataset size: {os.path.getsize(OUTPUT_ANNOTATED_H5AD)/(1024**2):.1f} MB')
print('\nREADY FOR PHASE III: Causal Discovery & Prior Triangulation!')

