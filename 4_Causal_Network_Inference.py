"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 4: Reconstruction of the triangulated tripartite causal Gene Regulatory Network (TF -> Enhancer -> Target Gene) integrating physical priors and interventional responses.

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 3]
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm

BASE_DIR = r'O:\Thesis\TNBC_Work\TNBC_Final_Dataset'
if not os.path.exists(BASE_DIR):
    BASE_DIR = os.getcwd()

PRIORS_FILE = os.path.join(BASE_DIR, 'TF_Enhancer_Priors.csv')
RESPONSES_FILE = os.path.join(BASE_DIR, 'Enhancer_TargetGene_Responses.csv')
OUTPUT_GRN_CSV = os.path.join(BASE_DIR, 'TNBC_Causal_Tripartite_GRN.csv')
FIG_DIR = os.path.join(BASE_DIR, 'manuscript_figures')
TAB_DIR = os.path.join(BASE_DIR, 'manuscript_tables')
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

print(f'Working directory: {BASE_DIR}')
print(f'Priors file: {PRIORS_FILE} (Exists: {os.path.exists(PRIORS_FILE)})')
print(f'Responses file: {RESPONSES_FILE} (Exists: {os.path.exists(RESPONSES_FILE)})')

# %% [Code Cell 5]
print('1. Loading ENCODE TF -> Enhancer priors...')
df_priors = pd.read_csv(PRIORS_FILE)
print(f'Loaded {len(df_priors):,} physical TF -> Enhancer binding events across {df_priors["tf_name"].nunique()} TFs.')

print('\n2. Loading Enhancer -> Target Gene perturbation responses...')
df_resp = pd.read_csv(RESPONSES_FILE)
print(f'Loaded {len(df_resp):,} functional response edges across {df_resp["enhancer_region"].nunique():,} perturbed enhancers.')

shared_enhancers = set(df_priors['enhancer_region']).intersection(set(df_resp['enhancer_region']))
print(f'\nShared regulatory enhancers between both layers: {len(shared_enhancers):,}')

# %% [Code Cell 7]
print('Triangulating physical TF binding with functional perturbation responses...')

# Subset to shared enhancers
priors_sub = df_priors[['tf_name', 'enhancer_region', 'confidence_score']].copy()
priors_sub.rename(columns={'confidence_score': 'chip_confidence'}, inplace=True)

resp_sub = df_resp[['enhancer_region', 'target_gene', 'log2_fc', 't_stat', 'effect_direction']].copy()

# Merge on enhancer_region
df_tripartite = pd.merge(priors_sub, resp_sub, on='enhancer_region', how='inner')
print(f'Constructed raw tripartite network with {len(df_tripartite):,} candidate paths.')

# Compute Causal Triangulation Score
# Score integrates physical binding strength (chip_confidence) with functional response magnitude
df_tripartite['causal_weight'] = np.abs(df_tripartite['log2_fc']) * np.abs(df_tripartite['t_stat']) * (df_tripartite['chip_confidence'] / 1000.0)
df_tripartite['causal_weight'] = df_tripartite['causal_weight'].round(4)

print('\nRaw Tripartite Path Sample:')
display(df_tripartite.head(4))

# %% [Code Cell 9]
print('Applying causal confidence thresholds...')
causal_mask = (np.abs(df_tripartite['t_stat']) >= 3.5) & (np.abs(df_tripartite['log2_fc']) >= 0.35)
df_causal_grn = df_tripartite[causal_mask].copy()
df_causal_grn.sort_values(by='causal_weight', ascending=False, inplace=True)

# Annotate evidence type
df_causal_grn['evidence_type'] = 'Physical_ChIP_seq_AND_Functional_CRISPRi'

print(f'Final High-Confidence Causal GRN contains: {len(df_causal_grn):,} directed edges.')
print(f'Regulating TFs: {df_causal_grn["tf_name"].nunique():,}')
print(f'Intermediate Enhancers: {df_causal_grn["enhancer_region"].nunique():,}')
print(f'Target Genes: {df_causal_grn["target_gene"].nunique():,}')

display(df_causal_grn.head(10))

# %% [Code Cell 11]
# 1. Save master Causal GRN table
df_causal_grn.to_csv(OUTPUT_GRN_CSV, index=False)
print(f'Saved Master Causal GRN to: {OUTPUT_GRN_CSV}')

# 2. Export Phase III Manuscript Table (Top Master TFs and their targets)
tf_summary = df_causal_grn.groupby('tf_name').agg(
    Controlled_Enhancers=('enhancer_region', 'nunique'),
    Downstream_Target_Genes=('target_gene', 'nunique'),
    Total_Causal_Edges=('target_gene', 'count'),
    Mean_Causal_Weight=('causal_weight', 'mean'),
    Max_Causal_Weight=('causal_weight', 'max')
).reset_index().sort_values(by='Total_Causal_Edges', ascending=False).head(20)
tf_summary['Mean_Causal_Weight'] = tf_summary['Mean_Causal_Weight'].round(2)
tf_summary['Max_Causal_Weight'] = tf_summary['Max_Causal_Weight'].round(2)

table3_xlsx = os.path.join(TAB_DIR, 'Table_3_Top_Causal_Master_Regulators.xlsx')
tf_summary.to_excel(table3_xlsx, index=False)
print(f'Saved Table 3 (Excel): {table3_xlsx}')

# 3. Export Publication Figure 4: Causal Network Architecture
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['figure.titlesize'] = 0
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))

# (A) Top 15 Master Causal TFs
top15_tfs = tf_summary.head(15)
sns.barplot(x=top15_tfs['Total_Causal_Edges'], y=top15_tfs['tf_name'], hue=top15_tfs['tf_name'], legend=False, ax=axes[0], palette='Reds_r', edgecolor='#333333', linewidth=0.6)
axes[0].set_xlabel('Number of Directed Causal Outgoing Edges')
axes[0].set_ylabel('Master Transcription Factor')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# (B) Distribution of Causal Weights
sns.histplot(df_causal_grn['causal_weight'], bins=40, ax=axes[1], color='#b2182b', edgecolor='#333333', linewidth=0.5)
axes[1].set_xlabel('Triangulated Causal Regulatory Weight')
axes[1].set_ylabel('Edge Frequency')
axes[1].set_yscale('log')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
fig5_png = os.path.join(FIG_DIR, 'Figure_4_Causal_Network_Architecture.png')
fig5_tiff = os.path.join(FIG_DIR, 'Figure_4_Causal_Network_Architecture.tiff')
fig.savefig(fig5_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(fig5_tiff, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved Figure 4: {fig5_png} & {fig5_tiff}')

print('\nPHASE III COMPLETE! Causal Directed Regulatory Network is fully constructed and triangulated.')

