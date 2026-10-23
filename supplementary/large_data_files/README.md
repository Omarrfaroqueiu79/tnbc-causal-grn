# Large-Scale High-Throughput Supplementary Datasets

This archive contains the complete, high-resolution functional genomic datasets supporting the manuscript:
**"Applying Causal Machine Learning and Multi-Omics Integration to Unravel Gene Regulatory Networks in Triple-Negative Breast Cancer (TNBC)"**

Public GitHub Repository: https://github.com/Omarrfaroqueiu79/tnbc-causal-grn

---

### Files Inventory:

1. **`Supplementary_Data_1_Dataset_and_QC_Summary.xlsx`** (12 KB)
   - Quality control metrics (median 6,519 UMIs, 2,138 genes, 12 Leiden clusters), cell filtering thresholds, and cohort parameters across 95,678 cells.
   - Location: `supplementary/Supplementary_Data_1_Dataset_and_QC_Summary.xlsx`

2. **`Supplementary_Table_1_Full_15k_TF_Enhancer_Priors.xlsx`** (710 KB)
   - Complete physical chromatin occupancy prior matrix: all 15,059 physical binding events across 161 unique transcription factors mapped to 819 targeted super-enhancers from ENCODE Factorbook ChIP-seq (`encRegTfbsClustered`).
   - Location: `supplementary/tables_excel/Supplementary_Table_1_Full_15k_TF_Enhancer_Priors.xlsx`

3. **`Supplementary_Data_2_Perturbation_Responses.csv.gz`** (25.1 MB)
   - Pairwise functional enhancer perturbation response matrix containing 1,110,566 candidate regulatory responses (|log2 FC| >= 0.25, |t| >= 2.5, p < 0.012) across 819 targeted enhancers and 32,799 target genes.
   - Location: `supplementary/large_data_files/Supplementary_Data_2_Perturbation_Responses.csv.gz`

4. **`Supplementary_Data_3` — Complete Reconstructed Tripartite Causal Network (6,366,350 paths)**
   - The full multi-million-edge causal regulatory network is partitioned into three self-contained compressed archives (<70 MB each for seamless GitHub hosting):
     - `Supplementary_Data_3_Part1_Causal_GRN.csv.gz` (52.1 MB, rows 1 to 2,122,116)
     - `Supplementary_Data_3_Part2_Causal_GRN.csv.gz` (62.5 MB, rows 2,122,117 to 4,244,232)
     - `Supplementary_Data_3_Part3_Causal_GRN.csv.gz` (65.7 MB, rows 4,244,233 to 6,366,350)
   - Also preserved as a unified monolithic file in local archives (`Supplementary_Data_3_Tripartite_Causal_GRN.csv.gz`, 180.2 MB).

---

### How to Load Supplementary Data 3 in Python:

```python
import pandas as pd

# Load all 3 parts into a single DataFrame:
df1 = pd.read_csv("Supplementary_Data_3_Part1_Causal_GRN.csv.gz")
df2 = pd.read_csv("Supplementary_Data_3_Part2_Causal_GRN.csv.gz")
df3 = pd.read_csv("Supplementary_Data_3_Part3_Causal_GRN.csv.gz")

full_grn = pd.concat([df1, df2, df3], ignore_index=True)
print("Reconstructed Tripartite Causal GRN shape:", full_grn.shape)
# Expected output: (6366350, 9)
```

### How to Merge in Bash / Linux:

```bash
gzip -dc Supplementary_Data_3_Part1_Causal_GRN.csv.gz > Supplementary_Data_3_Tripartite_Causal_GRN.csv
gzip -dc Supplementary_Data_3_Part2_Causal_GRN.csv.gz | tail -n +2 >> Supplementary_Data_3_Tripartite_Causal_GRN.csv
gzip -dc Supplementary_Data_3_Part3_Causal_GRN.csv.gz | tail -n +2 >> Supplementary_Data_3_Tripartite_Causal_GRN.csv
gzip -9 Supplementary_Data_3_Tripartite_Causal_GRN.csv
```
