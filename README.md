# Transcriptomic Analysis of Triple-Negative Breast Cancer

## Overview

This project presents a reproducible transcriptomic analysis of **Triple-Negative Breast Cancer (TNBC)** integrating differential expression, functional enrichment, candidate-gene prioritization, independent validation, molecular signature analysis, pathway-level analysis, ROC/AUC evaluation, protein-protein interaction analysis, and survival analysis.

The primary discovery analysis compares matched Normal breast tissue and TNBC samples from **GEO dataset GSE233242**.

The analysis identified two major transcriptional programs associated with the TNBC state:

1. **Increased cell-cycle and proliferative activity**
2. **Reduced lipid-associated and metabolic activity**

These programs were used to prioritize and evaluate a final ten-gene candidate panel:

**CDC20, BUB1, TRIP13, PLK1, AURKB, PNPLA2, PPARG, LIPE, LEP, and CIDEC**

The project is intended as a **computational and hypothesis-generating study**, rather than a clinical biomarker validation study.

---

## Research Question

> What transcriptional and biological pathways distinguish TNBC from Normal breast tissue, and which genes and molecular programs can be prioritized as candidate markers of these differences?

The analysis further evaluates whether the prioritized candidates show reproducible behavior in an independent cohort, whether they form coherent molecular networks, whether their expression can distinguish TNBC from Normal samples, and whether candidate-gene expression is associated with overall survival in a larger breast cancer cohort.

---

# Datasets

## Discovery Cohort

**GEO accession:** GSE233242

The discovery analysis used:

- 8 Normal samples
- 8 TNBC samples
- 16 samples total
- 39,376 genes/features

The analysis was restricted to matched Normal-TNBC samples.

## Independent Validation Cohort

**GEO accession:** GSE52194

The independent validation analysis contained:

- 3 Normal samples
- 5 TNBC samples

This cohort was used to evaluate whether the expression direction of the ten candidate genes and the two molecular programs observed in the discovery cohort could be reproduced externally.

## TCGA-BRCA Survival Cohort

Overall survival analysis was performed using the **TCGA-BRCA** cohort.

The survival analysis included:

- **1,073 patients** with both expression and survival data
- **150 observed deaths**
- **923 censored observations**

Expression data were obtained from the **UCSC Xena GDC TCGA-BRCA STAR-TPM** dataset.

---

# Analysis Workflow

```text
01 Data Exploration
        ↓
02 Differential Expression
        ↓
03 Functional Enrichment
        ↓
04 Candidate Prioritization
        ↓
05 Candidate Validation
        ↓
06 Visualization
        ↓
07 Independent Validation
        ↓
08 Molecular Signatures
        ↓
09 Hallmark GSEA
        ↓
10 Gene–Pathway Integration
        ↓
11 ROC/AUC Analysis
        ↓
12 Network Analysis
        ↓
13 Survival Analysis
        ↓
14 Drug–Target Analysis
        ↓
Integrated Interpretation
```

---

# Analysis Scripts

## 01. Data Exploration

`Scripts/01_explore_data.R`

- Loaded GEO sample metadata
- Identified Normal and TNBC samples
- Identified matched patients
- Loaded the raw RNA-seq count matrix
- Selected the final 16 samples
- Examined sequencing depth and count distributions

## 02. Differential Expression Analysis

`Scripts/02_differential_expression.R`

Differential expression analysis was performed using **DESeq2**.

The discovery analysis identified:

- **2,011 upregulated genes**
- **2,218 downregulated genes**

## 03. Functional Enrichment

`Scripts/03_functional_enrichment.R`

Functional enrichment was performed separately for upregulated and downregulated genes using:

- Gene Ontology Biological Process
- KEGG pathway analysis

The analysis identified strong enrichment of cell-cycle and mitotic processes among upregulated genes and lipid/metabolic pathways among downregulated genes.

## 04. Candidate Gene Prioritization

`Scripts/04_candidate_prioritization.R`

Candidate genes were prioritized using:

- Differential expression
- Statistical significance
- Pathway association
- Biological relevance

Ten final candidate genes were selected.

### Proliferation-associated genes

- CDC20
- BUB1
- TRIP13
- PLK1
- AURKB

### Lipid/metabolic-associated genes

- PNPLA2
- PPARG
- LIPE
- LEP
- CIDEC

## 05. Candidate Validation

`Scripts/05_candidate_validation.R`

The candidate genes were evaluated using:

- Expression-level comparison between Normal and TNBC
- Wilcoxon statistical testing
- Multiple-testing correction
- TNBC-specific Pearson correlation analysis

## 06. Visualization

`Scripts/06_visualization.R`

Generated visualizations including:

- PCA
- Differential-expression volcano plot
- GO enrichment plots
- KEGG enrichment plots
- Candidate-gene heatmap
- Candidate-gene boxplots
- TNBC candidate-gene correlation heatmap

## 07. Independent Validation

`Scripts/07_independent_validation.R`

The ten-gene panel was evaluated in the independent **GSE52194** cohort.

The original expression direction was reproduced for **8 of the 10 genes**.

### Direction reproduced

- CDC20
- BUB1
- TRIP13
- PLK1
- AURKB
- PNPLA2
- PPARG
- LIPE

### Direction not reproduced

- LEP
- CIDEC

The independent validation therefore provides **partial replication** of the candidate-gene expression program rather than complete replication of every candidate.

## 08. Molecular Signature Analysis

`Scripts/08_molecular_signatures.R`

The ten genes were organized into two predefined molecular programs.

| Molecular program | Genes |
|---|---|
| Proliferation | CDC20, BUB1, TRIP13, PLK1, AURKB |
| Metabolic | PNPLA2, PPARG, LIPE, LEP, CIDEC |

Signature scores were calculated for each sample.

### Discovery cohort

- Proliferation signature: **p = 0.0009391**
- Metabolic signature: **p = 0.0009391**

### Independent cohort

- Proliferation signature: **p = 0.03689**
- Metabolic signature: **p = 0.03689**

Both signatures therefore showed significant differences between Normal and TNBC samples in the analyzed cohorts.

## 09. Hallmark Gene Set Enrichment Analysis

`Scripts/09_GSEA.R`

Preranked Hallmark GSEA was performed using the DESeq2 test statistic as the ranking metric.

The analysis included:

- 24,750 ranked genes
- 50 Hallmark pathways
- `fgsea` for preranked enrichment analysis

Using an adjusted p-value threshold of 0.05:

- 12 pathways were enriched toward TNBC
- 14 pathways were enriched toward Normal tissue
- 26 pathways were significant overall

### Major TNBC-associated pathways

- G2M Checkpoint
- E2F Targets
- mTORC1 Signaling
- MYC Targets
- Glycolysis
- Mitotic Spindle

### Major Normal-associated pathways

- Adipogenesis
- Myogenesis
- Fatty Acid Metabolism
- Oxidative Phosphorylation
- Xenobiotic Metabolism
- Bile Acid Metabolism

These results provide pathway-level support for the proliferation-versus-metabolic interpretation of the candidate panel.

---

### 10. Gene–Pathway Integration

`Scripts/10_gene_pathway_integration.R`

The prioritized candidate genes were integrated with enriched biological pathways to examine how individual candidates relate to the major molecular programs identified in the analysis.

The analysis generated a candidate gene–pathway integration heatmap and related results linking the prioritized genes to enriched pathways.

---

# 11. ROC/AUC Analysis

`Scripts/11_ROC_AUC_Analysis.R`

ROC/AUC analysis was used to evaluate the ability of individual genes and molecular signatures to distinguish TNBC from Normal samples.

### Individual candidate genes

| Gene | AUC | Direction |
|---|---:|---|
| BUB1 | 0.942 | Higher in TNBC |
| CIDEC | 0.942 | Lower in TNBC |
| CDC20 | 0.937 | Higher in TNBC |
| LEP | 0.934 | Lower in TNBC |
| LIPE | 0.934 | Lower in TNBC |
| PNPLA2 | 0.915 | Lower in TNBC |
| PLK1 | 0.896 | Higher in TNBC |
| AURKB | 0.896 | Higher in TNBC |
| PPARG | 0.870 | Lower in TNBC |
| TRIP13 | 0.870 | Higher in TNBC |

The two molecular signatures showed:

- **Proliferation signature AUC = 1.00**
- **Metabolic signature AUC = 1.00**

The same AUC values were observed for both signatures in the independent GSE52194 cohort.

Because the discovery and independent cohorts are small, these AUC values should be interpreted as evidence of strong cohort separation rather than estimates of clinical diagnostic performance.

---

# 12. Protein-Protein Interaction and Network Analysis

`Scripts/12_network_analysis.R`

The ten candidate genes were evaluated using STRING-based protein-protein interaction information.

The high-confidence candidate network contained:

- **10 nodes**
- **18 edges**
- Average node degree = **3.6**
- Local clustering coefficient = **0.933**
- Expected number of edges = **1**
- PPI enrichment p-value = **3.13 × 10⁻¹⁴**

This indicates that the candidate genes form a highly interconnected network compared with the number of interactions expected by chance within the STRING framework.

### Network organization

The network contains two prominent molecular programs:

**Proliferation**

- CDC20
- BUB1
- TRIP13
- PLK1
- AURKB

**Metabolic**

- PNPLA2
- PPARG
- LIPE
- LEP
- CIDEC

Network analysis provides additional evidence that the selected genes are not an arbitrary collection of differentially expressed genes but participate in coherent interaction structures.

---

# 13. TCGA-BRCA Survival Analysis

`Scripts/13_Survival_Analysis.R`

The ten candidate genes were evaluated for association with overall survival using TCGA-BRCA.

The analysis included:

- **1,073 patients**
- **150 observed deaths**
- **923 censored observations**

For each gene, patients were divided into high- and low-expression groups using the median expression level.

The analysis included:

- Kaplan-Meier survival analysis
- Log-rank testing
- Cox proportional hazards regression
- 95% confidence intervals
- Benjamini-Hochberg FDR correction

## Main finding

PLK1 showed the strongest nominal association with overall survival:

- Hazard ratio = **1.46**
- 95% CI = **1.06–2.02**
- Cox p = **0.021**
- FDR = **0.214**

Thus, higher PLK1 expression was associated with increased mortality risk in the unadjusted Cox analysis, but this association did **not** remain statistically significant after correction for testing ten candidate genes.

AURKB and PNPLA2 showed suggestive trends but did not reach nominal statistical significance:

- AURKB: HR = 1.34, p = 0.073
- PNPLA2: HR = 0.75, p = 0.078

These survival findings should therefore be considered exploratory.

---

# 14. Drug–Target Analysis

`Scripts/14_Drug_Target_Analysis.R`

DGIdb-based drug–target analysis was performed to add a translational layer to the prioritized ten-gene panel.

The analysis:

- Queried the DGIdb GraphQL API for the candidate genes
- Retrieved reported drug–gene interactions
- Assigned candidate genes to the predefined proliferation or metabolic programs
- Removed incomplete and duplicate interaction records
- Summarized drug associations at the gene and drug levels
- Identified drugs associated with more than one prioritized candidate gene
- Generated candidate-gene drug-association and multi-target network visualizations

### Results

The analysis identified:

- **7 of the 10 candidate genes** with reported drug associations
- **472 unique drugs**
- **487 drug–gene interaction records**
- **7 multi-target drugs** associated with more than one prioritized candidate gene

### Gene-level drug associations

| Gene | Molecular program | Associated drugs | Interactions |
|---|---|---:|---:|
| PLK1 | Proliferation | 189 | 190 |
| PPARG | Metabolic | 169 | 173 |
| AURKB | Proliferation | 83 | 86 |
| LEP | Metabolic | 24 | 24 |
| BUB1 | Proliferation | 6 | 6 |
| LIPE | Metabolic | 5 | 5 |
| PNPLA2 | Metabolic | 3 | 3 |

PLK1 had the largest number of recorded drug associations, followed by PPARG and AURKB.

### Multi-target drugs

Seven compounds were associated with more than one prioritized candidate gene:

| Drug | Candidate genes | Molecular program(s) |
|---|---|---|
| AMORFRUTIN A | PLK1, PPARG | Proliferation, Metabolic |
| BENZBROMARONE | PLK1, PPARG | Proliferation, Metabolic |
| CHEMBL:CHEMBL184450 | PLK1, PPARG | Proliferation, Metabolic |
| ESTRADIOL VALERATE | LEP, PPARG | Metabolic |
| GW7647 | PLK1, PPARG | Proliferation, Metabolic |
| NVP-TAE684 | AURKB, PLK1 | Proliferation |
| TROGLITAZONE | LEP, PPARG | Metabolic |

The multi-target network highlights recurring relationships between the prioritized genes, particularly PLK1–PPARG and AURKB–PLK1.

These database-reported associations are **hypothesis-generating**. They do not establish therapeutic efficacy, drug sensitivity, or clinical suitability in TNBC.

### Drug–Target Outputs

- `Results/Drug_Target/`
- `Figures/Drug_Target/Candidate_Gene_Drug_Counts.png`
- `Figures/Drug_Target/Multi_Target_Drug_Network.png`

---

# Key Findings

## 1. Strong proliferative transcriptional program

TNBC samples showed increased expression of genes and pathways associated with:

- Mitotic division
- Chromosome segregation
- Cell-cycle checkpoints
- G2M checkpoint
- E2F targets
- MYC targets
- Mitotic spindle activity

The proliferation-associated candidates CDC20, BUB1, TRIP13, PLK1, and AURKB capture this program at the gene level.

## 2. Reduced lipid-associated and metabolic program

The TNBC samples showed reduced expression of genes and pathways associated with:

- Fatty acid metabolism
- PPAR signaling
- Lipolysis
- Adipogenesis
- Steroid metabolism
- Retinol metabolism
- Oxidative phosphorylation

The metabolic-associated candidates PNPLA2, PPARG, LIPE, LEP, and CIDEC capture this transcriptional program.

## 3. Independent replication is partial

The independent GSE52194 analysis reproduced the original expression direction for **8 of the 10 genes**.

LEP and CIDEC did not reproduce the discovery direction.

Therefore, the independent validation supports the overall molecular programs but does not provide complete gene-level replication.

## 4. Molecular signatures show strong cohort separation

Both the proliferation and metabolic signatures strongly separated TNBC from Normal samples in the discovery cohort and showed significant separation in the independent cohort.

This supports the use of the two programs as exploratory molecular summaries of the transcriptional differences observed in the analysis.

## 5. Candidate genes form a coherent interaction network

The ten genes formed a highly connected STRING network with 18 observed interactions compared with approximately 1 expected interaction.

The strong PPI enrichment supports functional coherence among the selected candidates.

## 6. Survival analysis identifies PLK1 as an exploratory prognostic candidate

PLK1 showed the strongest nominal association with overall survival in TCGA-BRCA.

However, its FDR-adjusted p-value was not significant.

Therefore:

> PLK1 should be considered an exploratory survival-associated candidate rather than a validated prognostic biomarker.

---

# Final Ten-Gene Candidate Panel

| Molecular program | Genes |
|---|---|
| Proliferation / mitosis | CDC20, BUB1, TRIP13, PLK1, AURKB |
| Lipid / metabolic | PNPLA2, PPARG, LIPE, LEP, CIDEC |

The panel captures two biologically contrasting transcriptional programs associated with the analyzed TNBC samples.

---

# Results Directory

The `Results/` directory contains:

- DESeq2 differential-expression results
- Upregulated and downregulated gene lists
- Annotated DEGs
- GO enrichment results
- KEGG enrichment results
- Candidate-gene prioritization results
- Final ten-gene shortlist
- Statistical validation results
- TNBC correlation matrix
- Independent validation results
- Molecular signature scores
- Hallmark GSEA results
- ROC/AUC results
- Network interaction results
- Network centrality results
- PPI enrichment results
- TCGA-BRCA survival results
- DGIdb drug–target interaction results
- Gene-level and drug-level DGIdb summaries
- Multi-target drug interaction results

---

# Figures Directory

The `Figures/` directory contains:

- PCA plots
- Volcano plots
- GO enrichment plots
- KEGG enrichment plots
- Candidate-gene heatmaps
- Candidate-gene boxplots
- Correlation heatmaps
- Molecular signature figures
- GSEA figures
- ROC/AUC figures
- PPI/network figures
- Kaplan-Meier survival plots
- Cox survival forest plot
- DGIdb candidate-gene drug-association plot
- Multi-target drug–candidate gene network

---

# Project Structure

```text
tnbc-transcriptomic-analysis/
│
├── Data/
│   ├── GSE233242_raw_counts_GRCh38.p13_NCBI.tsv.gz
│   ├── GSE233242_Sample_ID_to_GEO_ids.csv.gz
│   ├── Human.GRCh38.p13.annot.tsv.gz
│   ├── GSEA/
│   │   └── h.all.v2026.1.Hs.symbols.gmt
│   └── Independent_Validation/
│       └── E-GEOD-52194-raw-counts.tsv
│
├── Scripts/
│   ├── 01_explore_data.R
│   ├── 02_differential_expression.R
│   ├── 03_functional_enrichment.R
│   ├── 04_candidate_prioritization.R
│   ├── 05_candidate_validation.R
│   ├── 06_visualization.R
│   ├── 07_independent_validation.R
│   ├── 08_molecular_signatures.R
│   ├── 09_GSEA.R
│   ├── 10_gene_pathway_integration.R
│   ├── 11_ROC_AUC_Analysis.R
│   ├── 12_network_analysis.R
│   ├── 13_Survival_Analysis.R
│   └── 14_Drug_Target_Analysis.R
│
├── Results/
│   ├── Molecular_Signatures/
│   ├── GSEA/
│   ├── Independent_Validation/
│   ├── ROC_AUC/
│   ├── Network/
│   ├── Survival/
│   └── Drug_Target/
│
├── Figures/
│   ├── Molecular_Signatures/
│   ├── GSEA/
│   ├── Gene_Pathway_Integration/
│   ├── Network/
│   ├── Survival/
│   └── Drug_Target/
│
├── Interpretation/
│   └── biological_interpretation.md
│
├── README.md
├── .gitignore
└── Breast Cancer transcriptomics.Rproj
```

---

# Reproducibility

The analysis was performed in R.

Major packages used throughout the workflow include:

- DESeq2
- Biobase
- pheatmap
- clusterProfiler
- org.Hs.eg.db
- fgsea
- dplyr
- igraph
- survival
- survminer
- ggplot2
- httr2
- tidyr
- stringr

The drug–target analysis queries the **DGIdb GraphQL API** for reported drug–gene associations.

The discovery analysis is based on GSE233242.

The independent validation uses GSE52194.

Network analysis uses STRING interaction information.

Survival analysis uses TCGA-BRCA expression and survival data.

Large raw datasets are excluded from Git tracking where appropriate using `.gitignore`.

---

# Limitations

1. The discovery cohort contains only 16 samples: 8 Normal and 8 TNBC.
2. The independent validation cohort is also small, containing 3 Normal and 5 TNBC samples.
3. Independent validation reproduced the direction of 8/10 genes but not all candidates.
4. Molecular signature AUC values were calculated in small cohorts and should not be interpreted as estimates of clinical diagnostic performance.
5. ROC/AUC analyses describe discrimination within the analyzed datasets and do not establish clinical utility.
6. TNBC-specific correlation analysis was based on only eight TNBC samples and is exploratory.
7. Network analysis is based on STRING-derived interactions and does not establish physical interaction or causal relationships experimentally.
8. TCGA-BRCA survival analysis included 1,073 patients but only 150 observed deaths.
9. Survival analysis used median expression groups and was not adjusted for clinical covariates.
10. PLK1 showed a nominal survival association but did not remain significant after multiple-testing correction.
11. The ten genes should be considered candidate markers rather than clinically validated biomarkers.
12. The candidate panel should be considered hypothesis-generating rather than an established clinical molecular signature.
13. DGIdb associations represent database-reported drug–gene relationships and do not establish therapeutic efficacy or clinical suitability.
14. Only 7 of the 10 prioritized genes returned usable drug associations in the analyzed DGIdb query.
15. The number of recorded drug associations reflects database coverage and should not be interpreted as a measure of biological importance or therapeutic potential.

---

# Conclusion

This transcriptomic analysis identified extensive molecular differences between Normal and TNBC breast tissue.

Across differential expression, functional enrichment, GSEA, candidate prioritization, molecular signature analysis, independent validation, ROC/AUC analysis, network analysis, and survival analysis, the results consistently point toward two major transcriptional programs:

**Enhanced cell-cycle and proliferative activity**

and

**Reduced lipid-associated and metabolic activity**

The final ten-gene panel captures these contrasting programs:

**CDC20, BUB1, TRIP13, PLK1, AURKB, PNPLA2, PPARG, LIPE, LEP, and CIDEC.**

The independent validation reproduced the expression direction of 8 of the 10 genes and supported the broader molecular-program framework.

The two molecular signatures demonstrated strong separation between TNBC and Normal samples in both analyzed cohorts, while the candidate genes formed a highly interconnected STRING network.

In TCGA-BRCA survival analysis, PLK1 showed the strongest nominal association with overall survival. However, the association did not remain significant after multiple-testing correction.

DGIdb drug–target analysis added a translational layer, identifying reported drug associations for 7 candidate genes and highlighting PLK1, PPARG, and AURKB as candidates with substantial existing drug-association coverage. Several multi-target compounds connected candidate genes across the proliferation and metabolic programs. These findings are hypothesis-generating and do not establish therapeutic efficacy.

Together, these analyses provide a multi-level computational framework connecting differential expression, functional enrichment, candidate-gene prioritization, independent validation, molecular signatures, pathway integration, classification performance, network analysis, survival association, and drug–target relationships.

The drug–target analysis adds a translational hypothesis-generation layer, identifying substantial existing drug-association coverage for PLK1, PPARG, and AURKB and highlighting several multi-target compounds connecting the prioritized molecular programs.

The overall analysis remains exploratory and hypothesis-generating. Larger independent cohorts, experimental validation, and pharmacological studies are required before the candidate genes or associated compounds can be considered clinically actionable.

