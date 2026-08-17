# Transcriptomic Analysis of Triple-Negative Breast Cancer

## Overview

This project presents a transcriptomic analysis comparing Normal breast tissue with triple-negative breast cancer (TNBC) samples.

The analysis aims to identify differentially expressed genes, characterize the biological pathways associated with TNBC, and prioritize candidate genes that may represent distinct molecular programs associated with the disease.

The analysis identified two major transcriptional patterns:

- Increased cell-cycle and proliferative activity
- Reduced lipid-associated and metabolic activity

---

## Research Question

What transcriptional and biological pathways distinguish TNBC samples from Normal breast tissue, and which genes can be prioritized as candidate markers of these molecular differences?

---

## Dataset

The analysis used RNA-seq gene expression data consisting of:

- 8 Normal samples
- 8 TNBC samples
- 16 samples in the final analysis
- 39,376 genes/features

The analysis was restricted to matched Normal-TNBC samples.

---

## Analysis Workflow

The analysis was performed using R and was organized into six reproducible scripts.

### 1. Data Exploration and Sample Preparation

`Scripts/01_explore_data.R`

- Loaded GEO sample metadata
- Identified Normal and TNBC samples
- Identified matched patients
- Loaded the raw RNA-seq count matrix
- Selected the final 16 samples
- Examined sequencing depth and count distributions

### 2. Differential Expression Analysis

`Scripts/02_differential_expression.R`

Differential expression was performed using DESeq2.

The analysis identified:

- 2,011 upregulated genes
- 2,218 downregulated genes

### 3. Functional Enrichment

`Scripts/03_functional_enrichment.R`

Functional enrichment analysis was performed using:

- Gene Ontology (GO) Biological Process
- KEGG pathway analysis

The analysis was performed separately for upregulated and downregulated genes.

### 4. Candidate Gene Prioritization

`Scripts/04_candidate_prioritization.R`

Candidate genes were prioritized using differential expression, statistical significance, pathway association, and biological relevance.

Ten final candidates were selected.

### 5. Candidate Validation

`Scripts/05_candidate_validation.R`

The final candidates were evaluated using:

- Expression-level comparison between Normal and TNBC
- Wilcoxon statistical testing
- Multiple-testing correction
- TNBC-specific Pearson correlation analysis

### 6. Visualization

`Scripts/06_visualization.R`

The final analysis generated:

- PCA plot
- Volcano plot
- GO enrichment plots
- KEGG enrichment plots
- Candidate-gene heatmap
- Candidate-gene boxplots
- TNBC candidate-gene correlation heatmap

---

## Key Findings

### Differential Expression

The analysis identified extensive transcriptional differences between Normal and TNBC samples:

| Category | Number |
|---|---:|
| Upregulated genes | 2,011 |
| Downregulated genes | 2,218 |
| Total genes analyzed | 39,376 |

---

## Upregulated Molecular Program

The upregulated genes were strongly enriched for processes involving:

- Nuclear division
- Mitotic nuclear division
- Sister chromatid segregation
- Chromosome segregation
- Cell-cycle checkpoint signaling

The strongest GO Biological Process enrichment was nuclear division:

**Adjusted p-value = 1.63 × 10⁻²²**

The strongest KEGG pathway was:

**Cell cycle — adjusted p-value = 1.15 × 10⁻¹⁴**

These results indicate a strong proliferative and mitotic transcriptional program in the TNBC samples.

### Prioritized upregulated genes

- **CDC20**
- **BUB1**
- **TRIP13**
- **PLK1**
- **AURKB**

---

## Downregulated Molecular Program

Downregulated genes were enriched for metabolic, hormonal, and lipid-associated processes.

Important biological processes included:

- Fatty acid metabolic process
- Regulation of hormone levels
- Adaptive thermogenesis
- Retinol metabolism

KEGG analysis showed particularly strong enrichment of:

- **PPAR signaling — adjusted p-value = 2.58 × 10⁻⁵**
- **Regulation of lipolysis in adipocytes — adjusted p-value = 2.58 × 10⁻⁵**
- Steroid hormone biosynthesis
- Retinol metabolism
- Fatty acid degradation
- Adipocytokine signaling

These findings suggest altered lipid-associated and metabolic transcriptional programs in TNBC.

### Prioritized downregulated genes

- **PNPLA2**
- **PPARG**
- **LIPE**
- **LEP**
- **CIDEC**

---

## Final 10-Gene Candidate Panel

The final candidate panel consists of two major groups:

| Program | Genes |
|---|---|
| Proliferation / mitosis | CDC20, BUB1, TRIP13, PLK1, AURKB |
| Lipid / metabolic | PNPLA2, PPARG, LIPE, LEP, CIDEC |

The candidates showed significant expression differences between Normal and TNBC samples after multiple-testing correction.

---

## TNBC-Specific Correlation Analysis

Correlation analysis among the ten candidates was performed using the eight TNBC samples.

The proliferation-associated genes showed strong positive correlations, including:

- BUB1–TRIP13: r = 0.98
- TRIP13–PLK1: r = 0.97
- PLK1–AURKB: r = 0.95

The metabolic-associated genes also showed strong positive correlations:

- LIPE–CIDEC: r = 0.95
- LEP–CIDEC: r = 0.95
- PNPLA2–LIPE: r = 0.93

Several relationships between the proliferation and metabolic groups were negative.

This suggests that the ten candidates may represent two coordinated and opposing transcriptional programs within TNBC.

---

## Biological Interpretation

Taken together, the results suggest that the TNBC transcriptional state examined in this cohort is characterized by:

**Increased proliferative and mitotic activity**

alongside

**Reduced lipid-associated and metabolic activity.**

The coordinated expression of the prioritized genes provides additional support for this two-component molecular pattern.

The ten genes therefore represent a biologically coherent candidate molecular signature for further investigation.

A detailed interpretation of the results is provided in:

`Interpretation/biological_interpretation.md`

---

## Figures

### PCA

`Figures/PCA_Normal_vs_TNBC.png`

The PCA plot provides an overview of global transcriptomic variation between the Normal and TNBC samples.

### Differential Expression

`Figures/Volcano_Plot_TNBC_vs_Normal.png`

The volcano plot summarizes the magnitude and statistical significance of differential gene expression.

### Functional Enrichment

- `Figures/GO_BP_Upregulated.png`
- `Figures/GO_BP_Downregulated.png`
- `Figures/KEGG_Upregulated.png`
- `Figures/KEGG_Downregulated.png`

These figures summarize the biological processes and pathways associated with the upregulated and downregulated genes.

### Candidate Gene Expression

- `Figures/Final_10_Gene_Heatmap.png`
- `Figures/Final_10_Gene_Boxplots.png`

These figures show the expression patterns of the ten prioritized candidates.

### TNBC Gene Correlation

`Figures/TNBC_10_Gene_Correlation_Heatmap.png`

This figure shows gene-gene correlations among the ten candidates specifically within the TNBC samples.

---

## Results

The `Results/` directory contains:

- Complete DESeq2 results
- Upregulated and downregulated gene lists
- Annotated DEGs
- GO enrichment results
- KEGG enrichment results
- Candidate-gene prioritization results
- Final ten-gene shortlist
- Statistical validation results
- TNBC correlation matrix

---

## Limitations

The analysis has several limitations:

1. The final analysis included only 16 samples: 8 Normal and 8 TNBC.
2. The study identifies transcriptomic associations and does not establish causal relationships.
3. Candidate validation was performed within the same cohort rather than an independent dataset.
4. The TNBC correlation analysis was based on only 8 TNBC samples and should therefore be considered exploratory.
5. The ten genes should be considered candidate markers rather than clinically validated biomarkers.
6. Independent validation using larger cohorts and experimental studies is required.

---

## Reproducibility

The analysis is organized into sequential R scripts:

```text
01_explore_data.R
        ↓
02_differential_expression.R
        ↓
03_functional_enrichment.R
        ↓
04_candidate_prioritization.R
        ↓
05_candidate_validation.R
        ↓
06_visualization.R