# Supplementary Data Files Manifest

This directory contains the large-scale genomic matrices for the manuscript:
"Triangulating Single-Cell Perturbation Screens with Epigenetic Priors Reconstructs Causal Super-Enhancer Regulatory Networks in Triple-Negative Breast Cancer"

### Files Inventory:

1. **Supplementary_Data_1_Dataset_and_QC_Summary.xlsx**
   - Overview of single-cell sequencing parameters, QC cutoffs, and cohort definitions across 95,678 MDA-MB-231 cells.
   - Location: `Supplementary_Materials/Supplementary_Data_1_Dataset_and_QC_Summary.xlsx`

2. **Supplementary_Table_1_Full_15k_TF_Enhancer_Priors.xlsx**
   - 15,059 physical transcription factor binding events across 161 unique TFs mapped to targeted enhancers from ENCODE ChIP-seq.
   - Location: `Supplementary_Materials/Tables_Excel/Supplementary_Table_1_Full_15k_TF_Enhancer_Priors.xlsx`

3. **Supplementary_Data_2_Perturbation_Responses.csv.gz**
   - 1,110,566 candidate functional regulatory responses (|log2 FC| >= 0.25, |t| >= 2.5, p < 0.012) across 819 targeted super-enhancers and 32,799 genes.
   - Size: 25.09 MB (compressed GZIP CSV).
   - Location: `Supplementary_Materials/Large_Data_Files/Supplementary_Data_2_Perturbation_Responses.csv.gz`

4. **Supplementary_Data_3_Tripartite_Causal_GRN.csv.gz**
   - 6,366,350 triangulated high-confidence directed tripartite regulatory paths (TF -> Enhancer -> Target Gene) with composite causal regulatory weights.
   - Size: 180.24 MB (compressed GZIP CSV).
   - Location: `Supplementary_Materials/Large_Data_Files/Supplementary_Data_3_Tripartite_Causal_GRN.csv.gz`

5. **Supplementary_Display_1_Interactive_Network.html**
   - Standalone, interactive 2D/3D force-directed physics graph of the causal network.
   - Location: `Supplementary_Materials/Supplementary_Display_1_Interactive_Network.html`

All files are verified against the published manuscript text and equations.
