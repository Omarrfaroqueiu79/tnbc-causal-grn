# Point-by-Point Response to Supervisor Review Comments

**Manuscript Title:** *Integrating CRISPRi Perturbation Responses and Transcription-Factor Binding Priors for Regulatory Network Inference in Triple-Negative Breast Cancer*  
**Authors:** Omar Faroque¹, Mst. Homayra Jannaty¹, Shahriar Ahmed Abir², Suman Biswas¹\*  
*¹ Department of Statistics, Faculty of Sciences, Islamic University, Kushtia-7003, Bangladesh*  
*² Department of Biotechnology and Genetic Engineering, Faculty of Biological Sciences, Islamic University, Kushtia-7003, Bangladesh*  
*\* Corresponding Author:* `suman@stat.iu.ac.bd`

---

## Executive Overview of Manuscript Revisions

We express our sincere gratitude to the supervisor for this comprehensive, incisive, and constructive review. Every point raised has been systematically evaluated and integrated into the revised manuscript. All native mathematical equations (Equation 1, two-sample Welch $t$-test, PageRank formulation, XGBoost objective function, Cox proportional hazards regression) and formatting structures (Times New Roman, booktabs borders, pure black text) have been strictly preserved.

Below is the itemized summary of actions taken across the major themes:

---

### Theme 1: Title, Methodological Framing & Causal Language (Comments 1, 2, 43, 44, 48)

* **Comment 48 (Title Scoping):**  
  * *Action Taken:* We have fully adopted the supervisor's recommended title across all manuscript and supplementary files:  
    > **`Integrating CRISPRi Perturbation Responses and Transcription-Factor Binding Priors for Regulatory Network Inference in Triple-Negative Breast Cancer`**
* **Comments 1 & 2 (Causal Claims vs. Evidence Integration):**  
  * *Action Taken:* We revised the framing throughout the Abstract, Introduction, Results, and Discussion. The output is now explicitly described as a **"perturbation-supported, TF-binding-prior-informed regulatory network."** Direct causal claims are strictly reserved for the experimental enhancer perturbation links ($\text{Enhancer} \to \text{Gene}$ under CRISPRi), while the physical upstream link ($\text{TF} \to \text{Enhancer}$) is explicitly presented as an observational chromatin prior derived from ENCODE ChIP-seq.
* **Comments 43 & 44 (Tone & Terminology Moderation):**  
  * *Action Taken:* All promotional adjectives (*"catastrophic limitations,"* *"mathematically required,"* *"decisively overcomes,"* *"gold-standard"*) have been eliminated and replaced with objective, evidence-based academic language (*e.g., "substantial inferential challenges," "effectively reduces false discoveries," "established preclinical model"*).

---

### Theme 2: Network Decomposition & The 6.36M Paths (Comments 3, 32, 33)

* **Comments 3 & 32 (Network Entity Breakdown):**  
  * *Action Taken:* Section 3.4 and Table 3 have been updated to report the exact mathematical decomposition of the network:
    * **Total Tripartite Paths ($\text{TF} \to \text{Enhancer} \to \text{Gene}$):** **`6,366,349`**
    * **Unique (TF $\to$ Enhancer) Physical Edges:** **`5,935`**
    * **Unique (Enhancer $\to$ Gene) Functional Edges:** **`383,757`**
    * **Unique (TF $\to$ Target Gene) Functional Pairs:** **`542,760`**
    * **Unique Active Enhancers:** **`372`** (Mean targets/enhancer: $1,595.6$; Median: $1,596.0$)
    * **Unique Regulating TFs:** **`325`**
    * **Unique Downstream Target Genes:** **`2,504`**
* **Comment 33 (Direction of Effects & The "27 Downregulated Targets" in Table 4):**  
  * *Action Taken:* We investigated the underlying data and added a detailed technical footnote/explanation in Section 3.5. Under dCas9-KRAB super-enhancer repression, total cellular mRNA content decreases. Standard library depth normalization (CP10k) rescales the remaining counts, inducing a relative compositional upward shift (`pert_mean > ctrl_mean`), while only the 27 most invariant, highly expressed housekeeping and mitochondrial transcripts (`MALAT1`, `MT-CO2`, `BC_H1`–`BC_H9`) exhibit negative fold-changes. Because top master regulators bind across multiple shared active enhancers, this same invariant set emerges as the dominant downregulated subset across all leading drivers.

---

### Theme 3: Statistical Rigor & Multiple Testing (Comments 4, 5, 7, 16)

* **Comment 5 (Multiple-Testing FDR Correction):**  
  * *Action Taken:* Section 2.3 has been updated to explicitly state that the candidate perturbation tests were controlled using the **Benjamini–Hochberg False Discovery Rate (FDR)** procedure across the families of enhancer–gene hypotheses.
* **Comment 4 & 7 (Single-Cell Pseudoreplication & Guide Calling):**  
  * *Action Taken:* Section 2.1 now documents that single-cell guide calling mandated $\ge 5$ unique guide UMIs per cell. Aggregating cells per target enhancer locus (median ~14 sgRNAs and hundreds of cells per locus) effectively provided pseudobulk replicate averaging, mitigating droplet-level technical variance. The non-targeting control pool ($n > 5{,}000$ cells) is formally defined.
* **Comment 16 (Master Regulatory Score Formula):**  
  * *Action Taken:* Section 2.5 explicitly defines the normalization and weighting of the composite Master Regulatory Score:
    $$\text{MRS} = 0.5 \times \frac{\text{OutDegree}}{\max(\text{OutDegree})} + 0.5 \times \frac{\text{CumulativeWeight}}{\max(\text{CumulativeWeight})}$$

---

### Theme 4: Benchmarking Consistency & In Silico Validation (Comments 10, 11, 12, 14, 31)

* **Comments 11, 12 & 31 (Reconciliation of Discrepancies):**  
  * *Action Taken:* All performance metrics have been strictly harmonized across the Abstract, Results, Table 6, and Conclusion:
    * **Triangulated Framework:** **AUROC = 0.8999**, **AUPRC = 0.8745**, **Precision = 90.13%**, **FDR = 9.87%**.
    * **GENIE3:** **AUROC = 0.6698**, **AUPRC = 0.6323**, **Precision = 64.27%**, **FDR = 35.73%**.
    * **Observational Co-expression:** **AUROC = 0.7292**, **AUPRC = 0.6814**, **Precision = 67.73%**, **FDR = 32.27%**.
    * **Physical ChIP-seq Priors:** **AUROC = 0.8890**, **AUPRC = 0.7228**, **Precision = 75.60%**, **FDR = 24.40%**.
    * All unbenchmarked numbers (e.g., the previously mentioned 0.941/0.918) have been removed from the Conclusion.
* **Comment 14 (In Silico Knockout Simulation vs. TF Knockout):**  
  * *Action Taken:* Section 2.7 and Section 3.8 explicitly define the in silico knockout simulation as a **computational feature-silencing experiment ($do(\text{TF}=0)$)** within the learned regression model, avoiding any implication that wet-lab genetic TF knockout was performed.

---

### Theme 5: Clinical Survival & Therapeutic Framing (Comments 15, 19, 20, 40)

* **Comment 15 (Clinical Cohort Parameters):**  
  * *Action Taken:* Section 2.7 and Section 3.8 explicitly document the clinical cohort parameters:
    * **TCGA-BRCA ($n = 160$ TNBC patients):** 34 death events (21.3%), median follow-up 2.3 years.
    * **METABRIC ($n = 320$ TNBC patients):** 118 death events (36.9%), median follow-up 8.2 years.
    * Hazard ratios and 95% CIs reported for all 6 drivers (*JUND, MYC, FOS, EP300, MAX, CTCF*).
* **Comments 19, 20 & 40 (Therapeutic Network Disruption vs. Pharmacological Efficacy):**  
  * *Action Taken:* Language claiming "collapse of oncogenic core" or "chemo-refractory cure" was replaced with **"in silico topological network disruption"** (93% disruption under the triple epigenetic combination vs. 38% under chemotherapy expression signatures). Section 4.4 and Section 4.5 explicitly emphasize that these results represent computational hypotheses requiring future in vitro synergy assays (e.g., Chou–Talalay combination index).

---

### Theme 6: Expanded Limitations & Figure Evidential Clarity (Comments 28, 34, 41)

* **Comment 41 & 28 (Comprehensive Limitations Section):**  
  * *Action Taken:* Section 4.5 was expanded to systematically address:
    1. Preclinical boundaries of the MDA-MB-231 single cell-line model.
    2. Pan-cell-type ENCODE ChIP priors and the call for future breast-specific CUT&RUN profiles.
    3. Single-cell pseudoreplication and guide-assignment uncertainty.
    4. Model dependence of SHAP interpretations.
    5. The necessity of in vitro and in vivo pharmacological validation.
* **Comment 34 (Figure Captions):**  
  * *Action Taken:* Captions for Figures 1, 2, 7, and 8 were updated to provide complete evidential context (e.g., specifying that Figure 8 represents in silico topological network disruption and Figure 2 represents physical occupancy priors).

---

### Verification Summary of Updated Files
* **Main Manuscript Draft:** [`TNBC_Manuscript_Full_Draft.docx`](file:///L:/TNBC_Work/tnbc-causal-grn/manuscript/TNBC_Manuscript_Full_Draft.docx) & [`TNBC_Manuscript_Full_Draft.pdf`](file:///L:/TNBC_Work/tnbc-causal-grn/manuscript/TNBC_Manuscript_Full_Draft.pdf) (Clean 30 pages, zero blank pages).
* **Manuscript Tables:** [`Manuscript_Tables_All.docx`](file:///L:/TNBC_Work/Manuscript_Tables/Manuscript_Tables_All.docx) & [`Manuscript_Tables_All.pdf`](file:///L:/TNBC_Work/Manuscript_Tables/Manuscript_Tables_All.pdf) (Exactly 7 pages, Tables 1–7).
* **Supplementary Information:** [`Supplementary_Information.docx`](file:///L:/TNBC_Work/Supplementary_Materials/Supplementary_Information.docx) & [`Supplementary_Information.pdf`](file:///L:/TNBC_Work/Supplementary_Materials/Supplementary_Information.pdf).
* **Reporting Checklist:** [`STROBE_Statement_Checklist.docx`](file:///L:/TNBC_Work/Supplementary_Materials/STROBE_Statement_Checklist.docx) & [`STROBE_Statement_Checklist.pdf`](file:///L:/TNBC_Work/Supplementary_Materials/STROBE_Statement_Checklist.pdf).
