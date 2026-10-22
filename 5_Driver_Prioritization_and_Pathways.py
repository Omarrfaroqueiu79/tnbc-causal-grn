"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 5: Directed multigraph modeling, Google PageRank centrality, Master Regulatory Score (MRS) calculation, and MSigDB Hallmark pathway enrichment.

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 3]
import os
import numpy as np
import pandas as pd
import networkx as nx
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm

BASE_DIR = r'O:\Thesis\TNBC_Work\TNBC_Final_Dataset'
if not os.path.exists(BASE_DIR):
    BASE_DIR = os.getcwd()

GRN_FILE = os.path.join(BASE_DIR, 'TNBC_Causal_Tripartite_GRN.csv')
OUT_DRIVERS_CSV = os.path.join(BASE_DIR, 'Top_TNBC_Causal_Drivers.csv')
OUT_PATHWAYS_CSV = os.path.join(BASE_DIR, 'Enriched_TNBC_Phenotype_Pathways.csv')
FIG_DIR = os.path.join(BASE_DIR, 'manuscript_figures')
TAB_DIR = os.path.join(BASE_DIR, 'manuscript_tables')

print(f'Working directory: {BASE_DIR}')
print(f'Causal GRN input: {GRN_FILE} (Exists: {os.path.exists(GRN_FILE)})')

# %% [Code Cell 5]
print('Loading Causal GRN...')
# Reading relevant columns for network centrality
df_grn = pd.read_csv(GRN_FILE, usecols=['tf_name', 'target_gene', 'causal_weight', 'effect_direction'])
print(f'Loaded {len(df_grn):,} causal regulatory edges.')

# Aggregate to direct TF -> Target Gene regulatory summary
print('\nAggregating to direct TF -> Gene causal edges...')
tf_gene_edges = df_grn.groupby(['tf_name', 'target_gene']).agg(
    cumulative_causal_weight=('causal_weight', 'sum'),
    path_count=('causal_weight', 'count'),
    dominant_direction=('effect_direction', lambda x: x.mode()[0])
).reset_index()

print(f'Aggregated into {len(tf_gene_edges):,} unique TF -> Gene directed interactions.')
display(tf_gene_edges.head(5))

# %% [Code Cell 7]
print('Constructing NetworkX directed graph...')
G = nx.DiGraph()
for _, row in tqdm(tf_gene_edges.iterrows(), total=len(tf_gene_edges), desc='Adding Graph Edges'):
    G.add_edge(row['tf_name'], row['target_gene'], weight=row['cumulative_causal_weight'])

print(f'Graph Nodes: {G.number_of_nodes():,}, Edges: {G.number_of_edges():,}')

print('Computing PageRank scores...')
pagerank_scores = nx.pagerank(G, weight='weight')

# Compute summary metrics per TF
tf_list = list(df_grn['tf_name'].unique())
driver_metrics = []

for tf in tf_list:
    out_deg = G.out_degree(tf) if tf in G else 0
    pr = pagerank_scores.get(tf, 0.0)
    sub = tf_gene_edges[tf_gene_edges['tf_name'] == tf]
    cum_weight = sub['cumulative_causal_weight'].sum()
    
    driver_metrics.append({
        'tf_name': tf,
        'downstream_targets': out_deg,
        'pagerank_score': pr,
        'cumulative_causal_weight': round(cum_weight, 2),
        'upregulated_targets': (sub['dominant_direction'] == 'upregulated').sum(),
        'downregulated_targets': (sub['dominant_direction'] == 'downregulated').sum()
    })

df_drivers = pd.DataFrame(driver_metrics)
# Composite Master Regulatory Score (Normalized Out-degree * Cumulative Causal Weight)
norm_deg = df_drivers['downstream_targets'] / df_drivers['downstream_targets'].max()
norm_wt = df_drivers['cumulative_causal_weight'] / df_drivers['cumulative_causal_weight'].max()
df_drivers['master_regulatory_score'] = (0.5 * norm_deg + 0.5 * norm_wt).round(4)
df_drivers.sort_values(by='master_regulatory_score', ascending=False, inplace=True)
df_drivers.reset_index(drop=True, inplace=True)

print('\n--- Top 15 Prioritized Master Causal Regulators in TNBC ---')
display(df_drivers.head(15))

df_drivers.to_csv(OUT_DRIVERS_CSV, index=False)
print(f'Saved ranked drivers to: {OUT_DRIVERS_CSV}')

# %% [Code Cell 9]
print('Running MSigDB Hallmark pathway enrichment on top drivers...')
top_drivers = df_drivers.head(5)['tf_name'].tolist()
pathway_records = []

for tf in tqdm(top_drivers, desc='Enriching Top Drivers'):
    targets = tf_gene_edges[tf_gene_edges['tf_name'] == tf]['target_gene'].tolist()
    try:
        enr = gp.enrichr(
            gene_list=targets,
            gene_sets=['MSigDB_Hallmark_2020'],
            organism='human',
            outdir=None
        )
        res = enr.results
        res = res[res['Adjusted P-value'] < 0.05].copy()
        res['master_tf'] = tf
        pathway_records.append(res[['master_tf', 'Term', 'Adjusted P-value', 'Overlap', 'Combined Score', 'Genes']])
    except Exception as e:
        print(f'Skipped online enrichr for {tf}: {e}')

if pathway_records:
    df_pathways = pd.concat(pathway_records, ignore_index=True)
    df_pathways.rename(columns={'Term': 'Hallmark_Pathway', 'Adjusted P-value': 'FDR_Pval'}, inplace=True)
    df_pathways.to_csv(OUT_PATHWAYS_CSV, index=False)
    print(f'\nSaved enriched pathways to: {OUT_PATHWAYS_CSV}')
    display(df_pathways.head(10))
else:
    print('No online pathway results returned; check network connection.')

# %% [Code Cell 11]
# 1. Save Master Regulators Table to Excel
table4_xlsx = os.path.join(TAB_DIR, 'Table_4_Prioritized_TNBC_Drivers.xlsx')
df_drivers.head(25).to_excel(table4_xlsx, index=False)
print(f'Saved Table 4 (Excel): {table4_xlsx}')

# 2. Render Publication Figure: Driver Centrality & Hallmark Signatures
from adjustText import adjust_text

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['figure.titlesize'] = 0
plt.rcParams['axes.edgecolor'] = '#222222'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 9

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8))

# (a) Master Regulatory Scores of Top 15 Drivers
top15 = df_drivers.head(15)
sns.barplot(x=top15['master_regulatory_score'], y=top15['tf_name'], hue=top15['tf_name'], legend=False,
            ax=axes[0], palette='viridis', edgecolor='#222222', linewidth=0.6)
axes[0].set_xlabel('Master Regulatory Score')
axes[0].set_ylabel('Candidate Driver (TF)')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].text(-0.16, 1.06, 'a', transform=axes[0].transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# (b) Distribution of Target Genes Regulated per Driver (with collision avoidance)
sns.scatterplot(data=df_drivers, x='downstream_targets', y='cumulative_causal_weight',
                size='pagerank_score', sizes=(20, 200), color='#2b5c8f', alpha=0.65, ax=axes[1], legend=False)
axes[1].set_xlabel('Downstream Regulated Target Genes')
axes[1].set_ylabel('Cumulative Causal Regulatory Weight')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].text(-0.14, 1.06, 'b', transform=axes[1].transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left')

# Expand axes limits slightly to prevent clipping of labels
axes[1].set_xlim(df_drivers['downstream_targets'].min() - 50, df_drivers['downstream_targets'].max() + 160)
axes[1].set_ylim(-10000, df_drivers['cumulative_causal_weight'].max() * 1.12)

# Collect text annotations and auto-repel collisions using adjustText
texts = []
for _, row in df_drivers.head(6).iterrows():
    t = axes[1].text(
        row['downstream_targets'],
        row['cumulative_causal_weight'],
        row['tf_name'],
        fontsize=8,
        fontweight='bold',
        color='#111111'
    )
    texts.append(t)

adjust_text(
    texts,
    ax=axes[1],
    arrowprops=dict(arrowstyle='->', color='#444444', lw=0.6),
    avoid_self=True,
    force_points=0.3,
    force_text=0.5
)

plt.tight_layout()
fig_p4_png = os.path.join(FIG_DIR, 'Figure_5_Master_Driver_Prioritization.png')
fig_p4_tiff = os.path.join(FIG_DIR, 'Figure_5_Master_Driver_Prioritization.tiff')
fig.savefig(fig_p4_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(fig_p4_tiff, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Exported Figure 5: {fig_p4_png} & {fig_p4_tiff} (with anti-collision labels)')

print('\nPHASE IV COMPLETE! Key TNBC transcriptional drivers prioritized and mapped to phenotypic pathways.')

