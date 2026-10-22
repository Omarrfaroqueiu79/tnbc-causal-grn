"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 6: Non-linear TreeExplainer SHAP attribution modeling (XGBoost) and interactive force-directed graph generation via PyVis.

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 3]
import os
import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from sklearn.model_selection import train_test_split
import networkx as nx
from pyvis.network import Network
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

BASE_DIR = r'O:\Thesis\TNBC_Work\TNBC_Final_Dataset'
if not os.path.exists(BASE_DIR):
    BASE_DIR = os.getcwd()

GRN_FILE = os.path.join(BASE_DIR, 'TNBC_Causal_Tripartite_GRN.csv')
DRIVERS_FILE = os.path.join(BASE_DIR, 'Top_TNBC_Causal_Drivers.csv')
FIG_DIR = os.path.join(BASE_DIR, 'manuscript_figures')
TAB_DIR = os.path.join(BASE_DIR, 'manuscript_tables')
OUT_HTML = os.path.join(BASE_DIR, 'TNBC_Causal_Network_Interactive.html')
OUT_SHAP_CSV = os.path.join(TAB_DIR, 'Table_5_SHAP_Perturbation_Attribution.xlsx')

print(f'Working directory: {BASE_DIR}')
print(f'Causal GRN: {GRN_FILE} (Exists: {os.path.exists(GRN_FILE)})')
print(f'Drivers list: {DRIVERS_FILE} (Exists: {os.path.exists(DRIVERS_FILE)})')

# %% [Code Cell 5]
print('Loading Causal Network sample for Explainable AI...')
# Load representative causal paths across top drivers
df_grn = pd.read_csv(GRN_FILE, nrows=250000)
df_drivers = pd.read_csv(DRIVERS_FILE)

top_tfs = df_drivers.head(10)['tf_name'].tolist()
grn_top = df_grn[df_grn['tf_name'].isin(top_tfs)].copy()
print(f'Filtered to {len(grn_top):,} paths regulated by top 10 master drivers.')

# Pivot to construct TF -> Gene causal feature matrix
pivot_matrix = grn_top.pivot_table(index='target_gene', columns='tf_name', values='causal_weight', fill_value=0.0)

# Target variable: Target gene responsiveness (mean absolute log2 fold change)
gene_responsiveness = grn_top.groupby('target_gene')['log2_fc'].apply(lambda x: np.mean(np.abs(x)))
y = gene_responsiveness.loc[pivot_matrix.index].values
X = pivot_matrix.values
feature_names = list(pivot_matrix.columns)

print(f'Constructed feature matrix: {X.shape[0]} target genes x {X.shape[1]} master regulators.')

# Train Gradient Boosted Decision Trees (XGBoost)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = xgb.XGBRegressor(n_estimators=100, max_depth=4, learning_rate=0.05, random_state=42)
model.fit(X_train, y_train)

r2_score = model.score(X_test, y_test)
print(f'Trained XGBoost Causal Response Model | Test R2: {r2_score:.3f}')

# %% [Code Cell 7]
print('Computing SHAP values using TreeExplainer...')
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)

# Summary of feature attributions
mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
df_shap = pd.DataFrame({
    'Master_Regulator': feature_names,
    'Mean_Absolute_SHAP_Impact': mean_abs_shap,
    'Feature_Weight_Sum': X.sum(axis=0)
}).sort_values(by='Mean_Absolute_SHAP_Impact', ascending=False).reset_index(drop=True)

print('\n--- Master Regulators Ranked by SHAP Explanatory Power ---')
display(df_shap)

# Export Table 5
df_shap.to_excel(OUT_SHAP_CSV, index=False)
print(f'Saved SHAP Attribution Table: {OUT_SHAP_CSV}')

# %% [Code Cell 9]
print('Building interactive tripartite subnetwork for top TNBC drivers...')

# Select core subnetwork (Top 5 TFs, Top 15 Enhancers, Top 30 Target Genes)
core_tfs = df_drivers.head(5)['tf_name'].tolist()
sub_grn = df_grn[df_grn['tf_name'].isin(core_tfs)].sort_values(by='causal_weight', ascending=False).head(120)

net = Network(height='700px', width='100%', directed=True, bgcolor='#ffffff', font_color='#222222')
net.force_atlas_2based(gravity=-60, central_gravity=0.005, spring_length=120, spring_strength=0.08)

# Color schemes for tripartite layers
# TF: Coral / Red | Enhancer: Purple / Mauve | Target Gene: Teal / Blue
added_nodes = set()

for _, row in sub_grn.iterrows():
    tf = row['tf_name']
    enh = row['enhancer_region']
    gene = row['target_gene']
    weight = float(row['causal_weight'])
    direction = row['effect_direction']
    
    # Add TF node
    if tf not in added_nodes:
        net.add_node(tf, label=tf, color='#b2182b', size=25, title=f'Master TF: {tf}', shape='diamond')
        added_nodes.add(tf)
    
    # Add Enhancer node
    if enh not in added_nodes:
        short_enh = enh.split(':')[0] + ':' + enh.split(':')[1][:5] + '...'
        net.add_node(enh, label=short_enh, color='#7570b3', size=15, title=f'Enhancer Locus: {enh}', shape='dot')
        added_nodes.add(enh)
        
    # Add Target Gene node
    if gene not in added_nodes:
        net.add_node(gene, label=gene, color='#1b9e77', size=18, title=f'Target Gene: {gene}', shape='ellipse')
        added_nodes.add(gene)
    
    # Edges: TF -> Enhancer (Physical Binding)
    net.add_edge(tf, enh, color='#7570b3', width=1.5, title='ENCODE Physical Binding')
    
    # Edges: Enhancer -> Target Gene (Causal Perturbation)
    edge_col = '#d95f02' if direction == 'upregulated' else '#2b5c8f'
    net.add_edge(enh, gene, color=edge_col, width=max(1.0, min(weight/8.0, 4.5)),
                 title=f'Causal Impact: {direction} (weight: {weight:.2f})')

net.write_html(OUT_HTML)
print(f'Exported Interactive Network Graph to: {OUT_HTML}')

# Export Cytoscape tables for publication and external exploration
cyto_edges_path = os.path.join(TAB_DIR, 'Cytoscape_Network_Edges.csv')
sub_grn.to_csv(cyto_edges_path, index=False)
print(f'Exported Cytoscape Edge Table: {cyto_edges_path}')

# %% [Code Cell 11]
print('Rendering Figure 6: Model Explainability & Causal Attribution (Nature standard)...')

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['figure.titlesize'] = 0
plt.rcParams['axes.edgecolor'] = '#222222'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 9

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))

# (a) SHAP Feature Importance of Master Regulators
sns.barplot(data=df_shap, x='Mean_Absolute_SHAP_Impact', y='Master_Regulator', hue='Master_Regulator', legend=False,
            ax=axes[0], palette='flare', edgecolor='#222222', linewidth=0.6)
axes[0].set_xlabel('Mean |SHAP Value| (Causal Attribution)')
axes[0].set_ylabel('Master Transcriptional Driver')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].text(-0.16, 1.06, 'a', transform=axes[0].transAxes, fontsize=11, fontweight='bold')

# (b) SHAP Value Distribution vs Feature Weight (Top Driver: CTCF/POLR2A)
top_feature_idx = 0
top_feature_name = df_shap.iloc[0]['Master_Regulator']
axes[1].scatter(X[:, top_feature_idx], shap_values[:, top_feature_idx], color='#b2182b', s=8, alpha=0.5, rasterized=True)
axes[1].axhline(0, color='#666666', linestyle='--', linewidth=0.7)
axes[1].set_xlabel(f'Causal Regulatory Weight ({top_feature_name})')
axes[1].set_ylabel('SHAP Value (Impact on Regulation)')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].text(-0.14, 1.06, 'b', transform=axes[1].transAxes, fontsize=11, fontweight='bold')

plt.tight_layout()
fig6_png = os.path.join(FIG_DIR, 'Figure_6_Explainability_and_Causal_Pathways.png')
fig6_tiff = os.path.join(FIG_DIR, 'Figure_6_Explainability_and_Causal_Pathways.tiff')
fig.savefig(fig6_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(fig6_tiff, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Exported Figure 6: {fig6_png} & {fig6_tiff}')

print('\nPHASE V COMPLETE! Model explainability and interactive causal network visualization finished.')

