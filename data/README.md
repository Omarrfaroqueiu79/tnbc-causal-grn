# Data Acquisition and Preparation Guide

This directory contains metadata and automated instructions for obtaining the primary raw datasets utilized by the **TNBC Causal Tripartite GRN** discovery framework.

---

## 1. Primary Single-Cell CRISPRi Perturb-seq Dataset

* **Source Study:** Lewis et al. (2024), *"Combinatorial CRISPR screening reveals non-coding regulatory hubs in human breast cancer."*
* **NCBI GEO Accession:** [GSE224986](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224986)
* **Biological Model:** Triple-Negative Breast Cancer cell line `MDA-MB-231`.
* **Experimental Scale:**
  * 95,678 high-quality single cells post-filtering.
  * 11,975 sgRNAs targeting 819 somatic mutation-associated super-enhancers (selected from a 1,475-locus library; [GSE223812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223812)).
  * 32,799 reliably expressed genes.

### Download Command:
```bash
# Download raw transcript and guide matrices from GEO
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224986/suppl/GSE224986_MDAMB231_SM_trans_df.h5
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224986/suppl/GSE224986_MDAMB231_SM_sgRNA_df.h5
```

---

## 2. Epigenetic Multi-Omics Transcription Factor Binding Priors

* **Database:** ENCODE Factorbook & Clustered Transcription Factor ChIP-seq (`encRegTfbsClustered`).
* **Source:** UCSC Genome Browser & ENCODE REST API (ENCODE Project Consortium, 2012; Wang et al., 2013).
* **Coverage:** 340 transcription factors across 130 human cell lines.
* **Extraction:** Physical binding peaks overlapping targeted enhancer windows ($\pm 1,000$ bp flanking window) yielding 15,059 high-confidence physical TF–enhancer binding edges across 161 unique transcription factors.

---

## 3. Independent Clinical Validation Cohorts

* **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA):**
  * $n = 1,098$ patient transcriptomic profiles with linked clinical survival outcomes.
  * Stratified by receptor status: Triple-Negative Breast Cancer ($n = 160$ TNBC cases).
* **Molecular Taxonomy of Breast Cancer International Consortium (METABRIC):**
  * $n = 1,904$ patient tumor profiles with long-term overall survival and relapse follow-up.
  * Acquired via cBioPortal ([cbioportal.org](https://www.cbioportal.org)).
