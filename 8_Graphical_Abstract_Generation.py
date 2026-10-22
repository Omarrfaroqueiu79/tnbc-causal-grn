"""
TNBC Causal Tripartite GRN Discovery Pipeline
Stage 8: High-resolution publication-grade Graphical Abstract rendering and visual schematic generation.

Authors: Omar Faroque, Mst. Homayra Jannaty, Shahriar Ahmed Abir, Suman Biswas*
Department of Statistics, Islamic University, Kushtia-7003, Bangladesh
*Corresponding: suman@stat.iu.ac.bd
"""

# %% [Code Cell 2]
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image
from IPython.display import display

# Publication typography setup
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.bf'] = 'Arial:bold'

print("[INFO] Publication parameters and drawing libraries initialized successfully.")


# %% [Code Cell 4]
fig = plt.figure(figsize=(17.5, 11.5), dpi=300, facecolor='#FFFFFF')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 17.5)
ax.set_ylim(0, 11.5)
ax.axis('off')

# =============================================================
# 1. TOP-LEFT: CURATED DISCOVERY COHORT & HARMONIZATION
# =============================================================

# Ellipse Header: Multi-Omics Raw Sources
src_capsule = patches.FancyBboxPatch(
    (0.6, 10.0), 4.2, 1.1,
    boxstyle="round,pad=0.01,rounding_size=0.5",
    facecolor="#EFF6FF", edgecolor="#2563EB", linewidth=1.8
)
ax.add_patch(src_capsule)

# Two inner badges
b1 = patches.FancyBboxPatch((0.85, 10.55), 1.7, 0.42, boxstyle="round,pad=0.01,rounding_size=0.15",
                            facecolor="#FFFFFF", edgecolor="#3B82F6", linewidth=1.2)
ax.add_patch(b1)
ax.text(1.7, 10.76, "NCBI GEO\n(GSE224986)", fontsize=7.5, fontweight='bold', color='#1E40AF', ha='center', va='center')

b2 = patches.FancyBboxPatch((2.85, 10.55), 1.7, 0.42, boxstyle="round,pad=0.01,rounding_size=0.15",
                            facecolor="#FFFFFF", edgecolor="#3B82F6", linewidth=1.2)
ax.add_patch(b2)
ax.text(3.7, 10.76, "ENCODE\n(Factorbook)", fontsize=7.5, fontweight='bold', color='#1E40AF', ha='center', va='center')

# Cohort text inside capsule
ax.text(2.7, 10.32, "Harmonized Single-Cell Cohort", fontsize=8.5, fontweight='bold', color='#1E3A8A', ha='center')
ax.text(2.7, 10.14, "N = 95,678 Cells (CRISPRi + ChIP)", fontsize=7.8, fontweight='bold', color='#2563EB', ha='center')

# Title below capsule
ax.text(2.7, 9.72, "Curated Discovery Cohort", fontsize=11.5, fontweight='bold', color='#0F172A', ha='center')
ax.text(2.7, 9.50, "Triple-Negative Breast Cancer (MDA-MB-231)", fontsize=8.8, fontweight='bold', color='#B45309', ha='center')
ax.text(2.7, 9.32, "Pooled Super-Enhancer & Somatic Mutation Tiling", fontsize=7.6, color='#475569', ha='center')

# Arrow down to Harmonization Box
ax.annotate('', xy=(2.7, 8.85), xytext=(2.7, 9.20),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.35", color="#0F172A", lw=2.0))

# Harmonization Box (4 Tiers)
box_harm = patches.FancyBboxPatch(
    (0.5, 7.35), 4.4, 1.45,
    boxstyle="round,pad=0.01,rounding_size=0.06",
    facecolor="#FFFFFF", edgecolor="#0D9488", linewidth=1.5
)
ax.add_patch(box_harm)

tiers = [
    ("Tier 1: Single-Cell QC", "(>200 Genes, <15% Mito, 32,799 Genes)", "#0284C7", "#E0F2FE"),
    ("Tier 2: sgRNA Stratification", "(3,930 Single-Guide vs 7,767 Controls)", "#0D9488", "#CCFBF1"),
    ("Tier 3: Factorbook ChIP Priors", "(15,059 Physical TF-Enhancer Edges)", "#4F46E5", "#EEF2FF"),
    ("Tier 4: Interventional Responses", "(1,110,566 Functional Edges, |t| >= 2.5)", "#BE123C", "#FFE4E6")
]
for t_i, (t_title, t_desc, t_c, t_bg) in enumerate(tiers):
    ty = 8.55 - t_i * 0.33
    t_rect = patches.Rectangle((0.6, ty - 0.11), 4.2, 0.26, facecolor=t_bg, edgecolor='none')
    ax.add_patch(t_rect)
    ax.text(0.7, ty + 0.02, t_title, fontsize=7.5, fontweight='bold', color=t_c, va='center')
    ax.text(2.15, ty + 0.02, t_desc, fontsize=7.0, color='#1E293B', va='center')

ax.text(2.7, 7.15, "Computational Multi-Omics Harmonization", fontsize=9.8, fontweight='bold', color='#0F172A', ha='center')
ax.text(2.7, 6.96, "Orthogonal Physical & Functional Filtering", fontsize=7.6, color='#64748B', ha='center')

# Arrow down to Heatmap
ax.annotate('', xy=(2.7, 6.45), xytext=(2.7, 6.85),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.35", color="#0F172A", lw=2.0))

# Heatmap: Harmonized Interventional Matrix
# Headers
bar_green = patches.Rectangle((0.5, 6.22), 2.2, 0.22, facecolor='#16A34A', edgecolor='none')
ax.add_patch(bar_green)
ax.text(1.6, 6.33, "Single-Guide (Interventional)", fontsize=7.0, fontweight='bold', color='#FFFFFF', ha='center', va='center')

bar_red = patches.Rectangle((2.7, 6.22), 2.2, 0.22, facecolor='#DC2626', edgecolor='none')
ax.add_patch(bar_red)
ax.text(3.8, 6.33, "Negative Controls (Null)", fontsize=7.0, fontweight='bold', color='#FFFFFF', ha='center', va='center')

# Simulated heatmap blocks
np.random.seed(42)
hm_data = np.random.randn(8, 24)
hm_data[:, :12] += 1.3
hm_data[:, 12:] -= 0.9
hm_y0 = 4.35
hm_h = 1.82
hm_w = 4.4
for row_i in range(8):
    for col_i in range(24):
        val = hm_data[row_i, col_i]
        c = plt.cm.Reds(min(1.0, val / 2.5)) if val > 0 else plt.cm.Blues(min(1.0, -val / 2.5))
        rect = patches.Rectangle((0.5 + col_i * (hm_w / 24), hm_y0 + (7 - row_i) * (hm_h / 8)),
                                 (hm_w / 24), (hm_h / 8), facecolor=c, edgecolor='#FFFFFF', lw=0.4)
        ax.add_patch(rect)

row_labels = ["CTCF Occ.", "POLR2A Occ.", "EP300 Occ.", "JUND Occ.", "MYC Occ.", "EMT Resp.", "E2F Resp.", "Hypoxia"]
for r_i, r_lbl in enumerate(row_labels):
    ax.text(0.44, hm_y0 + (7 - r_i + 0.5) * (hm_h / 8), r_lbl, fontsize=6.8, fontweight='bold', color='#374151', ha='right', va='center')

ax.text(2.7, 4.15, "Harmonized Interventional Matrix", fontsize=9.8, fontweight='bold', color='#0F172A', ha='center')
ax.text(2.7, 3.96, "1,110,566 Functional Interventional Edges", fontsize=7.6, color='#64748B', ha='center')

# Branching connector line from Heatmap to bottom banner
ax.plot([4.9, 5.3, 5.3], [5.25, 5.25, 3.8], color='#B45309', lw=2.2)
ax.annotate('', xy=(5.3, 3.65), xytext=(5.3, 3.9),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.35", color="#B45309", lw=2.2))

# =============================================================
# 2. TOP RIGHT: CLINICAL VALIDATION (KM SURVIVAL)
# =============================================================
# Upward branch from Data Harmonization
ax.plot([4.9, 5.3, 5.3], [8.1, 8.1, 10.3], color='#1E293B', lw=2.0)
ax.annotate('', xy=(5.6, 10.3), xytext=(5.3, 10.3),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.35", color="#1E293B", lw=2.0))

ax.text(5.75, 10.65, "Clinical Validation", fontsize=12, fontweight='bold', color='#B45309')
ax.text(5.75, 10.40, "TCGA & METABRIC Cohorts (N = 3,002)", fontsize=8.5, fontweight='bold', color='#0284C7')
ax.text(5.75, 10.20, "In Silico Patient Survival Re-Analysis", fontsize=7.5, fontstyle='italic', color='#64748B')

# Mini dataset card
card_clin = patches.FancyBboxPatch((5.75, 8.8), 2.5, 1.25, boxstyle="round,pad=0.01,rounding_size=0.06",
                                  facecolor="#FFFFFF", edgecolor="#CBD5E1", linewidth=1.2)
ax.add_patch(card_clin)
badge_c = patches.Rectangle((5.95, 9.45), 2.1, 0.45, facecolor="#F0FDF4", edgecolor="#16A34A", lw=1.0)
ax.add_patch(badge_c)
ax.text(7.0, 9.75, "TCGA-BRCA / METABRIC", fontsize=7.2, fontweight='bold', color='#15803D', ha='center')
ax.text(7.0, 9.55, "Combined N = 3,002 Patients", fontsize=6.6, color='#15803D', ha='center')
ax.text(7.0, 9.20, "Multivariate Cox Proportional Hazards", fontsize=6.8, color='#374151', ha='center')
ax.text(7.0, 8.95, "Overall Survival Readout", fontsize=6.8, color='#64748B', ha='center')

# Arrow from dataset card to KM plot
ax.annotate('', xy=(8.75, 9.45), xytext=(8.35, 9.45),
            arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.3", color="#1E293B", lw=1.6))

# Step-line Kaplan-Meier Plot
km_x0 = 9.1
km_y0 = 8.65
km_w = 3.0
km_h = 1.9

# Axes
ax.plot([km_x0, km_x0 + km_w], [km_y0, km_y0], color='#1F2937', lw=1.2)
ax.plot([km_x0, km_x0], [km_y0, km_y0 + km_h], color='#1F2937', lw=1.2)

km_t = np.array([0, 10, 18, 25, 36, 48, 60])
km_s_low = np.array([1.0, 0.94, 0.88, 0.82, 0.76, 0.72, 0.70])
km_s_high = np.array([1.0, 0.80, 0.65, 0.52, 0.41, 0.35, 0.32])

def draw_step(times, survs, col, lbl):
    pts_x = [km_x0 + (times[0]/60)*km_w]
    pts_y = [km_y0 + survs[0]*km_h]
    for i in range(1, len(times)):
        pts_x.append(km_x0 + (times[i]/60)*km_w)
        pts_y.append(pts_y[-1])
        pts_x.append(km_x0 + (times[i]/60)*km_w)
        pts_y.append(km_y0 + survs[i]*km_h)
    ax.plot(pts_x, pts_y, color=col, lw=2.2, label=lbl)

draw_step(km_t, km_s_low, '#16A34A', 'Low Regulon')
draw_step(km_t, km_s_high, '#DC2626', 'High Regulon')

# KM Axes Ticks & Labels
ax.text(km_x0, km_y0 - 0.16, "0", fontsize=7.0, color='#64748B', ha='center')
ax.text(km_x0 + km_w/2, km_y0 - 0.16, "30", fontsize=7.0, color='#64748B', ha='center')
ax.text(km_x0 + km_w, km_y0 - 0.16, "60", fontsize=7.0, color='#64748B', ha='center')
ax.text(km_x0 + km_w/2, km_y0 - 0.34, "Time (Months)", fontsize=7.6, fontweight='bold', color='#374151', ha='center')

ax.text(km_x0 - 0.1, km_y0, "0.0", fontsize=7.0, color='#64748B', ha='right', va='center')
ax.text(km_x0 - 0.1, km_y0 + km_h/2, "0.5", fontsize=7.0, color='#64748B', ha='right', va='center')
ax.text(km_x0 - 0.1, km_y0 + km_h, "1.0", fontsize=7.0, color='#64748B', ha='right', va='center')
ax.text(km_x0 - 0.34, km_y0 + km_h/2, "Overall Survival", fontsize=7.6, fontweight='bold', color='#374151', rotation=90, ha='center', va='center')

# Legend for KM plot
ax.plot([km_x0 + 1.6, km_x0 + 1.9], [km_y0 + km_h - 0.25, km_y0 + km_h - 0.25], color='#16A34A', lw=2.0)
ax.text(km_x0 + 2.0, km_y0 + km_h - 0.25, "Low Risk", fontsize=6.8, color='#16A34A', fontweight='bold', va='center')
ax.plot([km_x0 + 1.6, km_x0 + 1.9], [km_y0 + km_h - 0.50, km_y0 + km_h - 0.50], color='#DC2626', lw=2.0)
ax.text(km_x0 + 2.0, km_y0 + km_h - 0.50, "High Risk", fontsize=6.8, color='#DC2626', fontweight='bold', va='center')

# Clinical Summary Text on Far Right
ax.text(12.6, 10.65, "Independent Survival Benefit", fontsize=11.5, fontweight='bold', color='#B45309')
ax.text(12.6, 10.38, "Independent Patient Cohorts (N = 3,002)", fontsize=8.8, fontweight='bold', color='#1E293B')
ax.text(12.6, 10.10, r"Master Hazard Ratio: HR = 1.88, 95% CI [1.44–2.46] (p < 0.001)", fontsize=8.8, fontweight='bold', color='#0284C7')
ax.text(12.6, 9.82, "• High expression of AP-1 (JUND/FOS) & MYC drives lethal progression", fontsize=7.6, color='#1E293B')
ax.text(12.6, 9.58, "• Validated independent of tumor stage, nodal status, and patient age", fontsize=7.6, color='#1E293B')
ax.text(12.6, 9.34, "• EP300 (HR = 1.45, p = 0.006) confirms master epigenetic dependency", fontsize=7.6, color='#1E293B')

# CR/PR grounding badge
recist_box = patches.FancyBboxPatch((12.6, 8.65), 1.1, 0.48, boxstyle="round,pad=0.01,rounding_size=0.06",
                                    facecolor="#FFFFFF", edgecolor="#16A34A", lw=1.2)
ax.add_patch(recist_box)
ax.text(13.15, 8.89, "HR <= 1.0\nProtected", fontsize=6.8, fontweight='bold', color='#16A34A', ha='center', va='center')

ax.text(13.9, 8.98, "Clinical Risk Stratification Grounding", fontsize=8.5, fontweight='bold', color='#B45309')
ax.text(13.9, 8.75, "• Stratifies High-Risk Recurrence vs Durable Chemosensitivity in TNBC", fontsize=7.4, color='#4B5563')

# =============================================================
# 3. MIDDLE RIGHT: CELLULAR GROUNDING & HALLMARK PHENOTYPING
# =============================================================
# Downward branch to Cellular Grounding
ax.plot([5.3, 5.3], [8.1, 6.75], color='#1E293B', lw=2.0)
ax.annotate('', xy=(5.6, 6.75), xytext=(5.3, 6.75),
            arrowprops=dict(arrowstyle="->,head_width=0.35,head_length=0.35", color="#1E293B", lw=2.0))

ax.text(5.75, 7.10, "Cellular Grounding", fontsize=12, fontweight='bold', color='#B45309')
ax.text(5.75, 6.85, "MDA-MB-231 Cohort (N = 95,678 Cells)", fontsize=8.5, fontweight='bold', color='#0284C7')
ax.text(5.75, 6.65, "Single-Cell Transcriptional Phenotyping", fontsize=7.5, fontstyle='italic', color='#64748B')

# Mini dataset card for single cells
card_sc = patches.FancyBboxPatch((5.75, 5.0), 2.5, 1.45, boxstyle="round,pad=0.01,rounding_size=0.06",
                                 facecolor="#FFFFFF", edgecolor="#CBD5E1", linewidth=1.2)
ax.add_patch(card_sc)
badge_sc = patches.Rectangle((5.95, 5.85), 2.1, 0.45, facecolor="#FAF5FF", edgecolor="#A855F7", lw=1.0)
ax.add_patch(badge_sc)
ax.text(7.0, 6.15, "NCBI GEO: GSE224986", fontsize=7.2, fontweight='bold', color='#7E22CE', ha='center')
ax.text(7.0, 5.95, "N = 95,678 Single Cells", fontsize=6.6, color='#7E22CE', ha='center')

# Mini UMAP dots inside card
np.random.seed(99)
u_x1 = np.random.normal(6.5, 0.18, 25)
u_y1 = np.random.normal(5.4, 0.12, 25)
u_x2 = np.random.normal(7.4, 0.18, 25)
u_y2 = np.random.normal(5.4, 0.12, 25)
ax.scatter(u_x1, u_y1, s=7, color='#0284C7', alpha=0.8)
ax.scatter(u_x2, u_y2, s=7, color='#BE123C', alpha=0.8)
ax.text(7.0, 5.12, "8 Leiden Cell States (UMAP)", fontsize=6.8, fontweight='bold', color='#475569', ha='center')

# Arrow from single cell card to bar chart
ax.annotate('', xy=(8.75, 5.8), xytext=(8.35, 5.8),
            arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.3", color="#1E293B", lw=1.6))

# Horizontal Hallmark Bar Chart
bar_x0 = 9.1
bar_y0 = 4.8
bar_w = 3.0
bar_h = 1.9

hallmark_bars = [
    ("EMT (Metastasis)", 0.95, "p = 1.2 x 10^-9", "#16A34A"),
    ("E2F Targets (Cycle)", 0.88, "p = 4.5 x 10^-8", "#0284C7"),
    ("G2M Checkpoint", 0.78, "p = 8.1 x 10^-7", "#D97706"),
    ("Hypoxia Signaling", 0.65, "p = 3.4 x 10^-5", "#DC2626")
]

for b_i, (b_name, b_val, b_pval, b_col) in enumerate(hallmark_bars):
    by = 6.45 - b_i * 0.44
    ax.text(bar_x0, by + 0.10, b_name, fontsize=7.6, fontweight='bold', color='#1E293B', va='center')
    b_len = b_val * 2.1
    rect = patches.Rectangle((bar_x0, by - 0.18), b_len, 0.22, facecolor=b_col, edgecolor='none')
    ax.add_patch(rect)
    ax.text(bar_x0 + b_len + 0.1, by - 0.07, b_pval, fontsize=7.2, fontweight='bold', color=b_col, va='center')

# Cellular Phenotyping Text on Far Right
ax.text(12.6, 7.10, "Single-Cell Phenotypic Grounding", fontsize=11.5, fontweight='bold', color='#B45309')
ax.text(12.6, 6.85, "High-Throughput scRNA-Seq Causal Enrichment", fontsize=8.8, fontweight='bold', color='#0284C7')
ax.text(12.6, 6.60, "Lewis et al. (MDA-MB-231, N = 95,678 Cells)", fontsize=8.0, color='#475569')
ax.text(12.6, 6.30, "• Epithelial-Mesenchymal Transition (EMT): p = 1.2 x 10^-9 (Master Hubs)", fontsize=7.6, color='#1E293B')
ax.text(12.6, 6.06, "• E2F Targets & G2M Progression: p = 4.5 x 10^-8 (MYC/MAX Obligate Control)", fontsize=7.6, color='#1E293B')
ax.text(12.6, 5.82, "• Hypoxia Core Adaptability: p = 3.4 x 10^-5 (EP300 / HDAC2 Scaffolding)", fontsize=7.6, color='#1E293B')
ax.text(12.6, 5.58, "• Mechanistically confirms super-enhancers dictate aggressive TNBC states", fontsize=7.6, fontstyle='italic', color='#475569')

# =============================================================
# 4. BOTTOM SECTION: THE CAUSAL ML FRAMEWORK (4 LARGE CARDS)
# =============================================================

# Wide Header Banner
banner = patches.FancyBboxPatch(
    (0.5, 3.25), 16.5, 0.55,
    boxstyle="round,pad=0.01,rounding_size=0.06",
    facecolor="#FFFBEB", edgecolor="#F59E0B", linewidth=1.8
)
ax.add_patch(banner)
ax.text(8.75, 3.52, "In Silico Causal Machine Learning Framework",
        fontsize=13.0, fontweight='bold', fontstyle='italic', color='#B45309', ha='center', va='center')

# 4 Modular Cards Under Banner
c_w = 3.9
c_gap = 0.3
c_xs = [0.5 + i * (c_w + c_gap) for i in range(4)]
card_y = 0.35
card_h = 2.75

for i in range(4):
    ax.annotate('', xy=(c_xs[i] + c_w/2, 3.12), xytext=(c_xs[i] + c_w/2, 3.25),
                arrowprops=dict(arrowstyle="->,head_width=0.3,head_length=0.3", color="#B45309", lw=1.5))

# -------------------------------------------------------------
# CARD 1: CAUSAL DISCOVERY (TRIPARTITE DAG)
# -------------------------------------------------------------
c1_x = c_xs[0]
card1 = patches.FancyBboxPatch((c1_x, card_y), c_w, card_h, boxstyle="round,pad=0.01,rounding_size=0.04",
                              facecolor="#FFFFFF", edgecolor="#CBD5E1", lw=1.2)
ax.add_patch(card1)
ax.text(c1_x + c_w/2, card_y + card_h - 0.24, "Causal Discovery (Tripartite DAG)", fontsize=9.2, fontweight='bold', color='#0F172A', ha='center')

# Node Layers: TF -> Enhancer -> Target Gene
tf_c = ["CTCF", "EP300", "JUND", "MYC"]
e_c = ["chr6..", "chr8..", "chr19.."]
g_c = ["ANXA2", "CALR", "CFL1"]

for idx, tf in enumerate(tf_c):
    ty = 2.45 - idx * 0.28
    ax.plot(c1_x + 0.7, ty, marker='o', markersize=7, color='#4338CA')
    ax.text(c1_x + 0.60, ty, tf, fontsize=6.8, fontweight='bold', color='#374151', ha='right', va='center')

for idx, e in enumerate(e_c):
    ey = 2.35 - idx * 0.34
    ax.plot(c1_x + 1.95, ey, marker='s', markersize=7, color='#0D9488')
    ax.text(c1_x + 1.95, ey - 0.14, e, fontsize=6.5, fontweight='bold', color='#0D9488', ha='center')

for idx, g in enumerate(g_c):
    gy = 2.35 - idx * 0.34
    ax.plot(c1_x + 3.2, gy, marker='^', markersize=7, color='#BE123C')
    ax.text(c1_x + 3.30, gy, g, fontsize=6.8, fontweight='bold', color='#374151', ha='left', va='center')

# Connect edges
for ti in range(4):
    for ei in [ti % 3, (ti + 1) % 3]:
        ax.plot([c1_x + 0.7, c1_x + 1.95], [2.45 - ti * 0.28, 2.35 - ei * 0.34], color='#C7D2FE', lw=1.0)
for ei in range(3):
    for gi in [ei, (ei + 1) % 3]:
        ax.plot([c1_x + 1.95, c1_x + 3.2], [2.35 - ei * 0.34, 2.35 - gi * 0.34], color='#FECDD3', lw=1.0)

# Triangulation Formula
ax.text(c1_x + c_w/2, 1.35, r"$W_{\mathrm{causal}} = C_{\mathrm{ChIP}} \times |\log_2\mathrm{FC}| \times |t|$",
        fontsize=7.8, fontweight='bold', color='#1E293B', ha='center')

ax.text(c1_x + c_w/2, 1.05, "• Dual-Evidence Triangulation\n• Pearl's Randomized Do-Operator\n• Suppresses Observational Confounders",
        fontsize=7.0, color='#374151', ha='center')

ax.text(c1_x + c_w/2, 0.52, "6,366,350 Directed Causal Paths", fontsize=8.0, fontweight='bold', color='#0284C7', ha='center')

# -------------------------------------------------------------
# CARD 2: QUANTITATIVE BENCHMARKING (ZERO FDR)
# -------------------------------------------------------------
c2_x = c_xs[1]
card2 = patches.FancyBboxPatch((c2_x, card_y), c_w, card_h, boxstyle="round,pad=0.01,rounding_size=0.04",
                              facecolor="#FFFFFF", edgecolor="#CBD5E1", lw=1.2)
ax.add_patch(card2)
ax.text(c2_x + c_w/2, card_y + card_h - 0.24, "Quantitative Benchmarking", fontsize=9.2, fontweight='bold', color='#0F172A', ha='center')

# Precision-Recall / ROC plot
roc_x = c2_x + 0.65
roc_y = 1.60
roc_w = 2.6
roc_h = 1.0

ax.plot([roc_x, roc_x + roc_w], [roc_y, roc_y], color='#475569', lw=1.0)
ax.plot([roc_x, roc_x], [roc_y, roc_y + roc_h], color='#475569', lw=1.0)

# Curves
ax.plot([roc_x, roc_x, roc_x + roc_w], [roc_y, roc_y + roc_h, roc_y + roc_h], color='#DC2626', lw=2.0, label='Causal ML (1.00)')
g_x = np.linspace(roc_x, roc_x + roc_w, 25)
g_y = roc_y + roc_h * (1 - np.exp(-4 * (g_x - roc_x) / roc_w))
ax.plot(g_x, g_y, color='#0284C7', linestyle='--', lw=1.5, label='GENIE3 (0.79)')
c_y = roc_y + roc_h * (g_x - roc_x) / roc_w
ax.plot(g_x, c_y, color='#D97706', linestyle=':', lw=1.5, label='Co-expr (0.49)')

# Labels
ax.text(roc_x, roc_y - 0.14, "0.0", fontsize=6.5, color='#64748B', ha='center')
ax.text(roc_x + roc_w, roc_y - 0.14, "1.0", fontsize=6.5, color='#64748B', ha='center')
ax.text(roc_x + roc_w/2, roc_y - 0.24, "Recall (Sensitivity)", fontsize=7.0, fontweight='bold', color='#374151', ha='center')
ax.text(roc_x - 0.12, roc_y + roc_h/2, "Precision", fontsize=7.0, fontweight='bold', color='#374151', rotation=90, ha='center', va='center')

ax.text(c2_x + c_w/2, 1.15, "• Held-Out Experimental CRISPRi Ground Truth\n• Triangulated Causal ML: AUPRC = 1.000\n• GENIE3 FDR = 17.0% | Co-expression FDR = 50.5%",
        fontsize=7.0, color='#374151', ha='center')

ax.text(c2_x + c_w/2, 0.52, "ZERO False Discovery Rate (FDR = 0.0%)", fontsize=8.0, fontweight='bold', color='#16A34A', ha='center')

# -------------------------------------------------------------
# CARD 3: EXPLAINABLE AI (SHAP TREEEXPLAINER)
# -------------------------------------------------------------
c3_x = c_xs[2]
card3 = patches.FancyBboxPatch((c3_x, card_y), c_w, card_h, boxstyle="round,pad=0.01,rounding_size=0.04",
                              facecolor="#FFFFFF", edgecolor="#CBD5E1", lw=1.2)
ax.add_patch(card3)
ax.text(c3_x + c_w/2, card_y + card_h - 0.24, "Explainable AI (SHAP Dynamics)", fontsize=9.2, fontweight='bold', color='#0F172A', ha='center')

# Mini SHAP switch curve plot
shap_x = c3_x + 0.65
shap_y = 1.60
shap_w = 2.6
shap_h = 1.0

ax.plot([shap_x, shap_x + shap_w], [shap_y, shap_y], color='#475569', lw=1.0)
ax.plot([shap_x, shap_x], [shap_y, shap_y + shap_h], color='#475569', lw=1.0)

np.random.seed(12)
sx = np.linspace(0, 20, 25)
sy_n = np.random.normal(0, 0.04, 25)
sy = np.where(sx < 8.5, 0.15 + sy_n, np.where(sx < 12.0, 0.85 + sy_n, 0.25 + sy_n))
sy = np.clip(sy, 0.05, 0.95)

screen_sx = shap_x + (sx / 20.0) * shap_w
screen_sy = shap_y + sy * shap_h

ax.plot(screen_sx, screen_sy, color='#BE123C', lw=1.8)
ax.scatter(screen_sx, screen_sy, color='#BE123C', s=12, alpha=0.9)

s_th = shap_x + (10.0 / 20.0) * shap_w
ax.plot([s_th, s_th], [shap_y, shap_y + shap_h], color='#0F172A', linestyle='--', lw=1.0)
ax.text(s_th, shap_y + shap_h + 0.06, "Switch Threshold (>10.0)", fontsize=6.5, fontweight='bold', color='#0F172A', ha='center')

ax.text(shap_x + shap_w/2, shap_y - 0.24, "JUND Causal Weight", fontsize=7.0, fontweight='bold', color='#374151', ha='center')
ax.text(shap_x - 0.12, shap_y + shap_h/2, "SHAP", fontsize=7.0, fontweight='bold', color='#374151', rotation=90, ha='center', va='center')

ax.text(c3_x + c_w/2, 1.15, "• XGBoost Regressor (Test R^2 = 0.742)\n• JUND (|SHAP| = 0.0296) Master Driver\n• Decodes Non-Monotonic Oncogene Activation",
        fontsize=7.0, color='#374151', ha='center')

ax.text(c3_x + c_w/2, 0.52, "AP-1 / MYC Non-Linear Drivers", fontsize=8.0, fontweight='bold', color='#BE123C', ha='center')

# -------------------------------------------------------------
# CARD 4: COUNTERFACTUAL INTERVENTION & THERAPY
# -------------------------------------------------------------
c4_x = c_xs[3]
card4 = patches.FancyBboxPatch((c4_x, card_y), c_w, card_h, boxstyle="round,pad=0.01,rounding_size=0.04",
                              facecolor="#FFFFFF", edgecolor="#CBD5E1", lw=1.2)
ax.add_patch(card4)
ax.text(c4_x + c_w/2, card_y + card_h - 0.24, "Counterfactual Intervention", fontsize=9.2, fontweight='bold', color='#0F172A', ha='center')

tx_x = c4_x + 0.4
tx_y0 = 1.65
tx_w = 3.1

tx_data = [
    ("Chemotherapy", 0.38, "#94A3B8"),
    ("BETi Monotherapy", 0.58, "#0284C7"),
    ("Dual (BETi+HATi)", 0.81, "#0D9488"),
    ("Triple Combo", 0.93, "#047857")
]

for idx, (tname, tpct, tcol) in enumerate(tx_data):
    ty = 2.45 - idx * 0.27
    ax.text(tx_x, ty + 0.06, tname, fontsize=7.0, fontweight='bold' if idx == 3 else 'normal', color='#1E293B', va='center')
    blen = tpct * (tx_w - 1.1)
    rect = patches.Rectangle((tx_x, ty - 0.10), blen, 0.14, facecolor=tcol, edgecolor='none')
    ax.add_patch(rect)
    ax.text(tx_x + blen + 0.08, ty - 0.03, f"{int(tpct*100)}%", fontsize=7.2, fontweight='bold', color=tcol, va='center')

ax.text(c4_x + c_w/2, 1.15, "• Pearl's Level-3 Counterfactual Intervention\n• Exp 1: do(BETi=1) -> 58% | Exp 2: Dual -> 81%\n• Exp 3: do(BETi+HATi+HDACi) -> 93% Shutdown",
        fontsize=7.0, color='#374151', ha='center')

ax.text(c4_x + c_w/2, 0.52, "93% Systemic Regulatory Collapse", fontsize=8.0, fontweight='bold', color='#047857', ha='center')

# Output paths
output_dir = r"l:\TNBC_Work\Overalll_Work\manuscript_figures"
png_path = os.path.join(output_dir, "Graphical_Abstract_LargeReadable.png")
tiff_path = os.path.join(output_dir, "Graphical_Abstract_LargeReadable.tiff")

# Overwrite primary
primary_png = os.path.join(output_dir, "Graphical_Abstract.png")
primary_tiff = os.path.join(output_dir, "Graphical_Abstract.tiff")

plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='#FFFFFF')
plt.savefig(tiff_path, dpi=300, bbox_inches='tight', facecolor='#FFFFFF')
plt.savefig(primary_png, dpi=300, bbox_inches='tight', facecolor='#FFFFFF')
plt.savefig(primary_tiff, dpi=300, bbox_inches='tight', facecolor='#FFFFFF')
plt.close()

# %% [Code Cell 6]
# Display inline verification preview
png_path = r"l:\TNBC_Work\Overalll_Work\manuscript_figures\Graphical_Abstract_LargeReadable.png"
if os.path.exists(png_path):
    img = Image.open(png_path)
    print(f"[VERIFIED] Dimensions: {img.size[0]} x {img.size[1]} pixels | Resolution: {img.info.get('dpi')} DPI | Mode: {img.mode}")
    display(img.resize((1200, 785)))
else:
    print("[ERROR] Output figure not found.")


