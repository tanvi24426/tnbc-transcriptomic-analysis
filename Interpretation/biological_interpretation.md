# Biological Interpretation

## 1. Overview

The transcriptomic analysis compared 8 Normal and 8 triple-negative breast cancer (TNBC) samples, resulting in a final analysis cohort of 16 samples. Differential expression analysis identified 2,011 upregulated genes and 2,218 downregulated genes in TNBC.

The analysis revealed two major and opposing transcriptional programs associated with TNBC:
  
  1. Increased cell-cycle and proliferative activity
2. Reduced lipid-associated and metabolic activity

These findings were supported by differential expression analysis, Gene Ontology (GO) enrichment, KEGG pathway analysis, candidate-gene prioritization, expression-level validation, and TNBC-specific gene-gene correlation analysis.

---
  
  ## 2. Global Transcriptomic Differences
  
  ### PCA
  
  The principal component analysis (PCA) was used to examine global transcriptional variation between samples.

The PCA plot showed separation between the Normal and TNBC samples, indicating that cancer status contributes substantially to the overall transcriptomic variation observed in the dataset.

This supports the differential expression analysis by demonstrating that the distinction between Normal and TNBC samples is reflected across the broader transcriptome rather than being driven by only a small number of individual genes.

**Figure:** `Figures/PCA_Normal_vs_TNBC.png`

---
  
  ## 3. Differential Expression
  
  Differential expression analysis identified:
  
  - 2,011 upregulated genes
- 2,218 downregulated genes

The volcano plot demonstrated extensive transcriptional remodeling in TNBC, with significant genes distributed in both the upregulated and downregulated directions.

This indicates that TNBC is associated with both activation and suppression of distinct biological programs rather than generalized changes in gene expression.

**Figure:** `Figures/Volcano_Plot_TNBC_vs_Normal.png`

---
  
  ## 4. Upregulated Biological Programs
  
  ### GO Biological Process Enrichment
  
  The upregulated genes were strongly enriched for processes related to cell division and mitosis.

The most significant terms included:
  
  | Biological process | Adjusted p-value |
  |---|---:|
  | Nuclear division | 1.63 × 10^-22 |
  | Organelle fission | 3.34 × 10^-20 |
  | Mitotic nuclear division | 5.25 × 10^-20 |
  | Sister chromatid segregation | 2.55 × 10^-18 |
  | Nuclear chromosome segregation | 2.98 × 10^-18 |
  | Mitotic sister chromatid segregation | 3.99 × 10^-18 |
  | Regulation of nuclear division | 2.90 × 10^-17 |
  | Cell-cycle checkpoint signaling | 1.24 × 10^-15 |
  
  The strong enrichment of these processes indicates that the upregulated TNBC transcriptome is dominated by genes involved in mitotic progression, chromosome segregation, and cell-cycle regulation.

The enrichment of cell-cycle checkpoint signaling is particularly relevant because checkpoint mechanisms are required for accurate progression through mitosis and maintenance of chromosome integrity.

**Figure:** `Figures/GO_BP_Upregulated.png`

### KEGG Enrichment

The strongest upregulated KEGG pathway was:
  
  - Cell cycle — adjusted p-value = 1.15 × 10^-14

Other significant pathways included:
  
  - DNA replication — adjusted p-value = 3.16 × 10^-3
- Viral carcinogenesis — adjusted p-value = 2.71 × 10^-3
- Steroid biosynthesis — adjusted p-value = 3.27 × 10^-3
- Cellular senescence — adjusted p-value = 1.60 × 10^-2
- p53 signaling — adjusted p-value = 4.29 × 10^-2

The strong enrichment of the cell-cycle and DNA-replication pathways independently supports the GO findings and indicates increased representation of proliferative and mitotic programs in TNBC.

Disease-associated pathways such as viral carcinogenesis should not be interpreted as evidence of viral infection. Their enrichment reflects the presence of genes shared across cancer-associated molecular mechanisms.

**Figure:** `Figures/KEGG_Upregulated.png`

---
  
  ## 5. Downregulated Biological Programs
  
  ### GO Biological Process Enrichment
  
  The downregulated genes showed enrichment for metabolic, hormonal, circulatory, and lipid-associated biological processes.

Important enriched terms included:
  
  | Biological process | Adjusted p-value |
  |---|---:|
  | Blood circulation | 5.90 × 10^-17 |
  | Olefinic compound metabolic process | 2.03 × 10^-14 |
  | Muscle system process | 1.07 × 10^-11 |
  | Regulation of hormone levels | 3.20 × 10^-11 |
  | Primary alcohol metabolic process | 2.33 × 10^-10 |
  | Fatty acid metabolic process | 1.45 × 10^-9 |
  | Adaptive thermogenesis | 2.20 × 10^-7 |
  | Retinol metabolic process | 2.95 × 10^-7 |
  
  The enrichment of fatty-acid metabolism, hormone regulation, adaptive thermogenesis, and retinol metabolism suggests that the downregulated transcriptome is associated with altered metabolic and lipid-related functions.

**Figure:** `Figures/GO_BP_Downregulated.png`

### KEGG Enrichment

The strongest downregulated pathways included:
  
  | Pathway | Adjusted p-value |
  |---|---:|
  | Tyrosine metabolism | 2.58 × 10^-5 |
  | Drug metabolism - cytochrome P450 | 2.58 × 10^-5 |
  | PPAR signaling pathway | 2.58 × 10^-5 |
  | Regulation of lipolysis in adipocytes | 2.58 × 10^-5 |
  | Neuroactive ligand-receptor interaction | 1.39 × 10^-4 |
  | Steroid hormone biosynthesis | 7.62 × 10^-4 |
  | Retinol metabolism | 1.26 × 10^-3 |
  | Fatty acid degradation | 1.17 × 10^-2 |
  | Adipocytokine signaling pathway | 2.68 × 10^-2 |
  | Arachidonic acid metabolism | 3.10 × 10^-2 |
  
  The particularly strong enrichment of PPAR signaling and regulation of lipolysis in adipocytes provides evidence for altered lipid-handling and metabolic signaling in TNBC.

**Figure:** `Figures/KEGG_Downregulated.png`

---
  
  ## 6. Candidate Gene Prioritization
  
  Ten genes were prioritized based on differential expression, statistical significance, pathway association, and biological relevance.

### Upregulated candidates

| Gene | log2 Fold Change | Adjusted p-value |
  |---|---:|---:|
  | CDC20 | +3.93 | 3.43 × 10^-18 |
  | BUB1 | +3.96 | 6.84 × 10^-16 |
  | TRIP13 | +2.84 | 5.01 × 10^-14 |
  | PLK1 | +2.67 | 1.07 × 10^-10 |
  | AURKB | +3.32 | 6.38 × 10^-10 |
  
  These genes represent the proliferative and mitotic component of the candidate panel.

### Downregulated candidates

| Gene | log2 Fold Change | Adjusted p-value |
  |---|---:|---:|
  | PNPLA2 | -2.95 | 7.01 × 10^-48 |
  | PPARG | -2.97 | 1.10 × 10^-32 |
  | LIPE | -4.73 | 2.42 × 10^-26 |
  | LEP | -5.92 | 8.56 × 10^-25 |
  | CIDEC | -6.65 | 1.04 × 10^-19 |
  
  These genes represent the lipid-associated and metabolic component of the candidate panel.

---
  
  ## 7. Candidate-Gene Expression Patterns
  
  The final ten-gene heatmap demonstrated a consistent difference in expression between Normal and TNBC samples.

The five proliferation-associated genes showed higher expression in TNBC, whereas the five metabolic-associated genes showed lower expression.

This indicates that the candidate genes collectively capture the major transcriptional differences identified during the broader differential expression analysis.

**Figure:** `Figures/Final_10_Gene_Heatmap.png`

The individual boxplots further demonstrated the expression differences between the two groups.

**Figure:** `Figures/Final_10_Gene_Boxplots.png`

---
  
  ## 8. Statistical Validation of Candidate Genes
  
  The ten prioritized candidates were evaluated using Wilcoxon rank-sum testing.

After multiple-testing correction, all ten genes remained significantly differentially expressed between Normal and TNBC samples.

Adjusted p-values ranged approximately from 0.00188 to 0.00388.

This supports the consistency of the expression differences within the analyzed cohort.

However, this represents within-cohort validation rather than independent external validation because the same cohort was used for candidate discovery and validation.

---
  
  ## 9. TNBC-Specific Gene-Gene Correlation
  
  Correlation analysis was performed using the 8 TNBC samples to examine relationships among the ten prioritized genes.

The proliferation-associated genes showed strong positive correlations:
  
  - BUB1 - TRIP13: r = 0.98
- TRIP13 - PLK1: r = 0.97
- PLK1 - AURKB: r = 0.95
- CDC20 - AURKB: r = 0.92

The metabolic-associated genes also showed strong positive correlations:
  
  - LIPE - CIDEC: r = 0.95
- LEP - CIDEC: r = 0.95
- PNPLA2 - LIPE: r = 0.93
- PNPLA2 - CIDEC: r = 0.90

Several relationships between the proliferation and metabolic groups were negative. Examples include:
  
  - PLK1 - LEP: r = -0.64
- CDC20 - LEP: r = -0.56
- PLK1 - PNPLA2: r = -0.57

These results suggest that the ten genes may form two coordinated expression programs within TNBC:
  
  ### Proliferation-associated module
  
  CDC20, BUB1, TRIP13, PLK1, AURKB

### Metabolic-associated module

PNPLA2, PPARG, LIPE, LEP, CIDEC

The generally negative relationships between these modules suggest an opposing transcriptional pattern between proliferative and metabolic programs.

**Figure:** `Figures/TNBC_10_Gene_Correlation_Heatmap.png`

Because the correlation analysis was based on only eight TNBC samples, these relationships should be considered exploratory and should not be interpreted as evidence of direct regulatory interactions.

---
  
  ## 10. Integrated Biological Interpretation
  
  Taken together, the results identify two major and opposing transcriptional programs associated with TNBC.

### Increased proliferative program

The upregulation of CDC20, BUB1, TRIP13, PLK1, and AURKB is supported by strong enrichment of mitotic nuclear division, chromosome segregation, cell-cycle checkpoint signaling, the cell-cycle pathway, and DNA replication.

This indicates that the TNBC samples have a strong transcriptional signature associated with cell proliferation and mitotic activity.

### Reduced metabolic program

The downregulation of PNPLA2, PPARG, LIPE, LEP, and CIDEC is supported by enrichment of fatty-acid metabolism, PPAR signaling, regulation of lipolysis, steroid hormone biosynthesis, retinol metabolism, and adipocytokine signaling.

This suggests that the TNBC transcriptional state is accompanied by reduced representation of specific lipid-associated and metabolic programs.

### Integrated model

The overall pattern can therefore be summarized as:
  
  **TNBC**
  
  → Increased cell-cycle and mitotic activity

**AND**
  
  → Reduced lipid-associated and metabolic activity

The strong within-group correlations further suggest that these changes may represent coordinated transcriptional programs rather than independent gene-level events.

---
  
  ## 11. Limitations
  
  Several limitations should be considered when interpreting these findings.

1. The analysis included only 16 samples: 8 Normal and 8 TNBC.
2. The study is based on transcriptomic associations and therefore does not establish causal relationships.
3. Candidate-gene validation was performed within the same cohort and was not an independent external validation.
4. The correlation analysis was based on only 8 TNBC samples and should therefore be considered exploratory.
5. The ten genes should be considered prioritized candidate genes rather than clinically validated biomarkers.
6. Larger independent cohorts and experimental validation are required to determine the reproducibility and biological significance of the identified signature.

---
  
  ## 12. Conclusion
  
  The transcriptomic analysis identified extensive molecular differences between Normal and TNBC samples, revealing two major and opposing biological programs.

The first was characterized by increased cell-cycle and mitotic activity, supported by strong enrichment of nuclear division, chromosome segregation, cell-cycle checkpoint signaling, the KEGG cell-cycle pathway, and DNA replication. The second was characterized by reduced lipid-associated and metabolic activity, supported by enrichment of fatty-acid metabolism, PPAR signaling, regulation of lipolysis, steroid hormone biosynthesis, retinol metabolism, and adipocytokine signaling.

Ten candidate genes were prioritized from these programs:
  
  **CDC20, BUB1, TRIP13, PLK1, AURKB, PNPLA2, PPARG, LIPE, LEP, and CIDEC.**
  
  The consistent expression differences and coordinated expression patterns of these genes support their prioritization as a biologically coherent candidate molecular signature.

Overall, the findings suggest that the TNBC transcriptional state examined in this cohort is characterized by a shift toward enhanced proliferative activity accompanied by reduced lipid-associated metabolic activity.

The identified ten-gene panel provides a starting point for further investigation, but independent validation in larger cohorts and experimental studies would be required before these genes could be considered robust biomarkers or therapeutic targets.