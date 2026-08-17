# Transcriptomic Analysis of Triple-Negative Breast Cancer

## Overview

This project presents a reproducible RNA-seq transcriptomic analysis comparing matched Normal breast tissue and triple-negative breast cancer (TNBC) samples.

The analysis aims to identify differentially expressed genes, characterize biological pathways associated with TNBC, and prioritize candidate genes representing distinct molecular programs.

The analysis identified two major transcriptional patterns:

- **Increased cell-cycle and proliferative activity**
- **Reduced lipid-associated and metabolic activity**

---

## Research Question

> What transcriptional and biological pathways distinguish TNBC samples from Normal breast tissue, and which genes can be prioritized as candidate markers of these molecular differences?

---

## Dataset

The analysis used RNA-seq gene expression data from the **NCBI Gene Expression Omnibus (GEO) dataset GSE233242**.

The final analysis included:

- **8 Normal samples**
- **8 TNBC samples**
- **16 samples in total**
- **39,376 genes/features**

The analysis was restricted to matched Normal-TNBC samples.

---

# Analysis Workflow

The analysis was performed using **R** and organized into six sequential scripts.

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