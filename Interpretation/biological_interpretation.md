# Biological Interpretation

## 1. Overview

This project investigates transcriptomic differences between Normal breast tissue and triple-negative breast cancer (TNBC), integrating differential expression, functional enrichment, candidate-gene prioritization, independent validation, molecular signatures, Hallmark GSEA, ROC/AUC analysis, protein-protein interaction analysis, and survival analysis.

The discovery analysis compared 8 Normal and 8 TNBC samples, resulting in a final analysis cohort of 16 samples. Differential expression analysis identified 2,011 upregulated genes and 2,218 downregulated genes. fileciteturn14file0L5-L12

Across the analysis, two major and opposing transcriptional programs emerged:

1. **Increased cell-cycle and proliferative activity**
2. **Reduced lipid-associated and metabolic activity**

The final candidate panel consisted of:

**CDC20, BUB1, TRIP13, PLK1, AURKB, PNPLA2, PPARG, LIPE, LEP, and CIDEC.**

The biological interpretation below integrates the discovery cohort with the independent validation cohort, molecular signature analysis, pathway-level analysis, network analysis, and TCGA-BRCA survival analysis.

---

# 2. Global Transcriptomic Differences

## PCA

Principal component analysis (PCA) was used to examine global transcriptional variation between samples.

The PCA showed separation between the Normal and TNBC samples, indicating that cancer status contributes substantially to the overall transcriptomic variation observed in the discovery dataset.

This supports the differential expression analysis by showing that the distinction between Normal and TNBC is reflected across the broader transcriptome rather than being restricted to a small number of individual genes.

**Figure:** `Figures/PCA_Normal_vs_TNBC.png`

---

# 3. Differential Expression

Differential expression analysis identified:

- **2,011 upregulated genes**
- **2,218 downregulated genes**

The volcano plot demonstrated extensive transcriptional remodeling in TNBC, with significant genes distributed in both directions.

This indicates that the TNBC samples are characterized by coordinated activation and suppression of distinct biological programs rather than a generalized increase or decrease in gene expression.

**Figure:** `Figures/Volcano_Plot_TNBC_vs_Normal.png`

---

# 4. Upregulated Biological Programs

## GO Biological Process Enrichment

The upregulated genes were strongly enriched for processes related to cell division and mitosis.

Important enriched terms included:

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

The strong enrichment of these processes indicates that the TNBC transcriptome is strongly associated with mitotic progression, chromosome segregation, and cell-cycle regulation.

The enrichment of cell-cycle checkpoint signaling is particularly relevant because checkpoint mechanisms coordinate progression through the cell cycle and help maintain chromosome integrity.

**Figure:** `Figures/GO_BP_Upregulated.png`

## KEGG Enrichment

The strongest upregulated KEGG pathway was:

- **Cell cycle — adjusted p-value = 1.15 × 10^-14**

Other significant pathways included:

- DNA replication — adjusted p-value = 3.16 × 10^-3
- Viral carcinogenesis — adjusted p-value = 2.71 × 10^-3
- Steroid biosynthesis — adjusted p-value = 3.27 × 10^-3
- Cellular senescence — adjusted p-value = 1.60 × 10^-2
- p53 signaling — adjusted p-value = 4.29 × 10^-2

The strong enrichment of cell-cycle and DNA-replication pathways independently supports the GO findings and indicates increased representation of proliferative and mitotic programs in TNBC.

The enrichment of cancer-associated pathway categories such as viral carcinogenesis should not be interpreted as evidence of viral infection. Such pathway annotations can contain genes that participate in broader cancer-associated molecular mechanisms.

**Figure:** `Figures/KEGG_Upregulated.png`

---

# 5. Downregulated Biological Programs

## GO Biological Process Enrichment

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

The enrichment of fatty-acid metabolism, hormone regulation, adaptive thermogenesis, and retinol metabolism indicates that the downregulated transcriptome is associated with altered lipid-related and metabolic functions.

**Figure:** `Figures/GO_BP_Downregulated.png`

## KEGG Enrichment

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

The strong enrichment of PPAR signaling and regulation of lipolysis provides pathway-level support for altered lipid-handling and metabolic signaling in the TNBC samples.

**Figure:** `Figures/KEGG_Downregulated.png`

---

# 6. Candidate Gene Prioritization

Ten genes were prioritized based on differential expression, statistical significance, pathway association, and biological relevance.

## Upregulated candidates

| Gene | log2 Fold Change | Adjusted p-value |
|---|---:|---:|
| CDC20 | +3.93 | 3.43 × 10^-18 |
| BUB1 | +3.96 | 6.84 × 10^-16 |
| TRIP13 | +2.84 | 5.01 × 10^-14 |
| PLK1 | +2.67 | 1.07 × 10^-10 |
| AURKB | +3.32 | 6.38 × 10^-10 |

These genes represent the proliferative and mitotic component of the candidate panel.

## Downregulated candidates

| Gene | log2 Fold Change | Adjusted p-value |
|---|---:|---:|
| PNPLA2 | -2.95 | 7.01 × 10^-48 |
| PPARG | -2.97 | 1.10 × 10^-32 |
| LIPE | -4.73 | 2.42 × 10^-26 |
| LEP | -5.92 | 8.56 × 10^-25 |
| CIDEC | -6.65 | 1.04 × 10^-19 |

These genes represent the lipid-associated and metabolic component of the candidate panel.

---

# 7. Candidate-Gene Expression Patterns

The final ten-gene heatmap demonstrated a clear difference in expression between Normal and TNBC samples.

The five proliferation-associated genes showed higher expression in TNBC, whereas the five metabolic-associated genes showed lower expression.

This indicates that the candidate genes collectively capture the major transcriptional contrast identified during differential expression analysis.

**Figure:** `Figures/Final_10_Gene_Heatmap.png`

The individual boxplots further demonstrated the expression differences between the two groups.

**Figure:** `Figures/Final_10_Gene_Boxplots.png`

---

# 8. Statistical Validation of Candidate Genes

The ten prioritized candidates were evaluated using Wilcoxon rank-sum testing.

After multiple-testing correction, all ten genes remained significantly differentially expressed between Normal and TNBC samples, with adjusted p-values approximately ranging from 0.00188 to 0.00388.

This supports the consistency of the expression differences within the discovery cohort.

However, this analysis is **within-cohort validation**, because the same cohort was used for candidate discovery and statistical validation. It should therefore not be interpreted as independent external validation.

---

# 9. Independent Validation

The candidate panel was evaluated in the independent **GSE52194** cohort containing 3 Normal and 5 TNBC samples.

The original expression direction was reproduced for **8 of the 10 genes**:

- CDC20
- BUB1
- TRIP13
- PLK1
- AURKB
- PNPLA2
- PPARG
- LIPE

The direction was not reproduced for:

- LEP
- CIDEC

This provides partial independent replication of the candidate-gene program.

The result is important because it indicates that the overall transcriptional pattern is not completely restricted to the discovery dataset, while also demonstrating that individual candidates may behave differently across cohorts.

Because the independent cohort is very small, lack of replication for individual genes should be interpreted cautiously.

---

# 10. TNBC-Specific Gene-Gene Correlation

Correlation analysis was performed using the 8 TNBC samples to examine coordinated expression among the ten prioritized genes.

## Proliferation-associated correlations

Several proliferation-associated genes showed strong positive correlations:

- BUB1 - TRIP13: r = 0.98
- TRIP13 - PLK1: r = 0.97
- PLK1 - AURKB: r = 0.95
- CDC20 - AURKB: r = 0.92

## Metabolic-associated correlations

The metabolic-associated genes also showed strong positive correlations:

- LIPE - CIDEC: r = 0.95
- LEP - CIDEC: r = 0.95
- PNPLA2 - LIPE: r = 0.93
- PNPLA2 - CIDEC: r = 0.90

## Between-program relationships

Several relationships between the two groups were negative:

- PLK1 - LEP: r = -0.64
- CDC20 - LEP: r = -0.56
- PLK1 - PNPLA2: r = -0.57

These relationships are consistent with the interpretation that the candidate genes form two coordinated expression programs within the TNBC samples:

### Proliferation-associated module

**CDC20, BUB1, TRIP13, PLK1, AURKB**

### Metabolic-associated module

**PNPLA2, PPARG, LIPE, LEP, CIDEC**

However, because these correlations are based on only eight TNBC samples, they should be considered exploratory and should not be interpreted as evidence of direct regulatory interactions.

**Figure:** `Figures/TNBC_10_Gene_Correlation_Heatmap.png`

---

# 11. Molecular Signature Analysis

The ten genes were organized into two predefined molecular signatures:

| Molecular signature | Genes |
|---|---|
| Proliferation | CDC20, BUB1, TRIP13, PLK1, AURKB |
| Metabolic | PNPLA2, PPARG, LIPE, LEP, CIDEC |

Signature scores were calculated to summarize the activity of each molecular program rather than relying only on individual genes.

## Discovery cohort

Both signatures showed significant differences between Normal and TNBC samples:

- Proliferation signature: **p = 0.0009391**
- Metabolic signature: **p = 0.0009391**

## Independent cohort

The same two programs also showed significant separation:

- Proliferation signature: **p = 0.03689**
- Metabolic signature: **p = 0.03689**

The fact that both signatures remain significant in the independent cohort provides additional support for the broader molecular-program interpretation, despite the small sample size.

---

# 12. Hallmark GSEA

Preranked Hallmark GSEA was performed using the differential-expression test statistic as the ranking metric.

The analysis included:

- 24,750 ranked genes
- 50 Hallmark pathways
- `fgsea` for preranked enrichment analysis

Using an adjusted p-value threshold of 0.05:

- 12 pathways were enriched toward TNBC
- 14 pathways were enriched toward Normal tissue
- 26 pathways were significant overall

## Major TNBC-associated pathways

The strongest TNBC-associated Hallmark programs included:

- G2M Checkpoint
- E2F Targets
- mTORC1 Signaling
- MYC Targets
- Glycolysis
- Mitotic Spindle

These pathways provide an independent pathway-level view of the proliferative state identified by the candidate genes and GO/KEGG analyses.

## Major Normal-associated pathways

Major pathways enriched toward Normal tissue included:

- Adipogenesis
- Myogenesis
- Fatty Acid Metabolism
- Oxidative Phosphorylation
- Xenobiotic Metabolism
- Bile Acid Metabolism

These results provide additional support for the metabolic component of the interpretation.

The convergence between candidate-gene analysis, GO/KEGG enrichment, and Hallmark GSEA strengthens the interpretation that the observed differences represent coordinated biological programs rather than isolated gene-level changes.

---

# 13. ROC/AUC Analysis

ROC/AUC analysis was used to assess how well the candidate genes and molecular signatures separated TNBC from Normal samples.

## Individual candidate genes

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

The individual candidates therefore showed strong separation between the Normal and TNBC samples within the analyzed datasets.

The two molecular signatures showed:

- **Proliferation signature AUC = 1.00**
- **Metabolic signature AUC = 1.00**

The same AUC values were observed in the independent GSE52194 cohort.

These results indicate strong cohort-level discrimination. However, the very small sample sizes mean that the AUC values should **not** be interpreted as evidence of clinical diagnostic performance. Larger independent datasets would be required to estimate generalizable predictive performance.

---

# 14. Protein-Protein Interaction and Network Analysis

The ten candidate genes were evaluated using STRING-based protein-protein interaction information.

The high-confidence candidate network contained:

- **10 nodes**
- **18 edges**
- Average node degree = **3.6**
- Local clustering coefficient = **0.933**
- Expected number of edges = **1**
- PPI enrichment p-value = **3.13 × 10^-14**

The observed network therefore contained substantially more interactions than expected by chance within the STRING framework.

This supports functional coherence among the candidate genes.

The network can be broadly interpreted in terms of the same two biological programs:

### Proliferation network

**CDC20, BUB1, TRIP13, PLK1, AURKB**

These genes are connected to cell-cycle progression, mitotic regulation, chromosome segregation, and checkpoint-related processes.

### Metabolic network

**PNPLA2, PPARG, LIPE, LEP, CIDEC**

These genes represent the lipid-associated/metabolic component of the candidate panel.

The network analysis strengthens the argument that the candidate panel is biologically structured rather than simply a collection of statistically significant genes.

Importantly, STRING-based network connectivity does not demonstrate direct physical interaction or causal regulation.

**Figure:** `Figures/Network/Candidate_PPI_Network.png`

---

# 15. TCGA-BRCA Survival Analysis

The ten candidate genes were evaluated for association with overall survival using TCGA-BRCA.

The final matched survival-expression cohort included:

- **1,073 patients**
- **150 observed deaths**
- **923 censored observations**
- No missing survival time
- No missing event status

Patients were divided into high- and low-expression groups using median gene expression.

The analysis included:

- Kaplan-Meier survival analysis
- Log-rank testing
- Cox proportional hazards regression
- 95% confidence intervals
- Benjamini-Hochberg FDR correction

## Main finding: PLK1

PLK1 showed the strongest nominal association with overall survival:

- Hazard ratio = **1.46**
- 95% CI = **1.06–2.02**
- Cox p = **0.021**
- Cox FDR = **0.214**

This indicates that higher PLK1 expression was associated with increased mortality risk in the unadjusted Cox model.

However, the association did **not** remain statistically significant after multiple-testing correction across the ten candidate genes.

Therefore, PLK1 should be interpreted as an **exploratory survival-associated candidate**, rather than a validated prognostic biomarker.

## Other suggestive associations

AURKB showed:

- HR = 1.34
- Cox p = 0.073

PNPLA2 showed:

- HR = 0.75
- Cox p = 0.078

Neither reached nominal statistical significance.

The remaining candidates did not show strong evidence of association with overall survival in this analysis.

### Interpretation

The survival analysis adds a distinct layer to the project because it asks whether the candidate genes are associated not only with TNBC-versus-Normal transcriptional differences but also with outcome variation across the broader TCGA-BRCA cohort.

The nominal PLK1 association is biologically consistent with the broader proliferative program identified in the discovery analysis, but the lack of FDR significance means that this observation should remain hypothesis-generating.

The survival analysis was not adjusted for clinical covariates and should therefore not be interpreted as demonstrating independent prognostic value.

**Figures:** `Figures/Survival/`  
**Results:** `Results/Survival/10_Gene_Survival_Analysis.csv`

---


# 16. Drug–Target Analysis

A DGIdb-based drug–target analysis was performed to add a translational layer to the prioritized ten-gene panel.

The analysis identified reported drug associations for **7 of the 10 candidate genes**, resulting in **487 drug–gene interaction records involving 472 unique drugs**.

## Gene-level drug associations

The greatest number of reported drug associations was observed for:

| Gene | Molecular program | Number of associated drugs | Number of interactions |
|---|---|---:|---:|
| PLK1 | Proliferation | 189 | 190 |
| PPARG | Metabolic | 169 | 173 |
| AURKB | Proliferation | 83 | 86 |
| LEP | Metabolic | 24 | 24 |
| BUB1 | Proliferation | 6 | 6 |
| LIPE | Metabolic | 5 | 5 |
| PNPLA2 | Metabolic | 3 | 3 |

PLK1 therefore had the largest number of recorded drug associations, followed by PPARG and AURKB. This provides additional translational context for candidates that were already supported by the transcriptomic, network, and survival analyses.

## Multi-target drug associations

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

The multi-target network highlights two recurring relationships. **PLK1–PPARG** is represented by several compounds, connecting the proliferation and metabolic candidate programs, while **AURKB–PLK1** is represented by NVP-TAE684 within the proliferation program. LEP and PPARG are jointly represented by estradiol valerate and troglitazone.

These relationships do not imply that the compounds are effective treatments for TNBC. Rather, they identify existing database-reported drug–gene associations that can be used to generate hypotheses for future drug-repurposing or mechanistic studies.

## Translational interpretation

The drug–target analysis adds a translational layer to the candidate-gene framework. **PLK1 and PPARG** were the two genes with the largest numbers of recorded drug associations, while **AURKB** also showed substantial drug-association coverage.

PLK1 is particularly notable because it connects several independent analysis layers in this project: it is strongly upregulated in the discovery cohort, belongs to the proliferation-associated molecular program, participates in the candidate PPI network, showed strong cohort-level discrimination in ROC analysis, and showed the strongest nominal survival association in TCGA-BRCA. However, its survival association did not remain significant after FDR correction. The drug-association results therefore provide additional context rather than independent evidence of therapeutic value.

Similarly, PPARG showed extensive drug-association coverage and belongs to the metabolic program supported by downregulation, PPAR/lipolysis pathway enrichment, molecular signatures, and network analysis. The database associations provide a rationale for further investigation of PPARG-related pharmacological mechanisms but do not establish whether targeting PPARG would be beneficial in TNBC.

**Figures:**  
- `Figures/Drug_Target/Candidate_Gene_Drug_Counts.png`
- `Figures/Drug_Target/Multi_Target_Drug_Network.png`

**Results:** `Results/Drug_Target/`

# 17. Integrated Biological Model

Taken together, the analyses converge on two major and opposing transcriptional programs associated with the TNBC samples.

## Program 1: Increased proliferation and mitotic activity

The proliferation-associated genes:

**CDC20, BUB1, TRIP13, PLK1, AURKB**

are supported by:

- Strong differential upregulation
- GO enrichment for nuclear division and chromosome segregation
- KEGG enrichment for cell cycle and DNA replication
- Hallmark enrichment for G2M checkpoint and E2F targets
- Strong positive correlations among proliferation-associated genes
- High classification performance in ROC analysis
- A coherent PPI network

This provides multi-level evidence for an enhanced proliferative transcriptional state.

## Program 2: Reduced lipid-associated and metabolic activity

The metabolic-associated genes:

**PNPLA2, PPARG, LIPE, LEP, CIDEC**

are supported by:

- Strong differential downregulation
- GO enrichment for fatty-acid and metabolic processes
- KEGG enrichment for PPAR signaling and regulation of lipolysis
- Hallmark enrichment for adipogenesis and fatty-acid metabolism in the Normal direction
- Strong positive correlations among metabolic-associated genes
- High classification performance in ROC analysis
- Coherent network connectivity

This supports an altered lipid-associated metabolic transcriptional state in the TNBC samples.

## Integrated interpretation

The overall pattern can therefore be summarized as:

```text
                         TNBC
                          │
            ┌─────────────┴─────────────┐
            ↓                           ↓
   Increased proliferative       Reduced lipid-associated
        transcriptional              / metabolic program
           activity
            │                           │
     CDC20, BUB1, TRIP13,       PNPLA2, PPARG, LIPE,
       PLK1, AURKB                 LEP, CIDEC
            │                           │
       Cell cycle /                  PPAR / lipid
       mitosis / G2M               metabolism / adipogenesis
```

The two programs are not proposed here as direct causally opposing pathways. Rather, they represent two coordinated transcriptional patterns that distinguish the TNBC samples from the Normal samples in the analyzed datasets.

---

# 18. Evidence Across Analysis Layers

The interpretation is supported at multiple levels:

| Analysis layer | Main observation | Interpretation |
|---|---|---|
| PCA | Normal and TNBC samples separate | Global transcriptomic differences |
| DESeq2 | 2,011 upregulated and 2,218 downregulated genes | Extensive transcriptional remodeling |
| GO | Mitotic processes up; metabolic processes down | Functional separation |
| KEGG | Cell cycle up; PPAR/lipolysis down | Pathway-level support |
| Candidate genes | 5 proliferation + 5 metabolic genes | Focused molecular panel |
| Within-cohort validation | All 10 remain significant after correction | Consistent discovery-cohort signal |
| Independent validation | 8/10 genes reproduce direction | Partial external replication |
| Signatures | Both programs significant in both cohorts | Program-level replication |
| Hallmark GSEA | G2M/E2F/MYC vs adipogenesis/fatty-acid metabolism | Independent pathway support |
| ROC/AUC | Strong individual-gene discrimination | Strong cohort separation |
| PPI network | 18 edges vs 1 expected | Functional/network coherence |
| Survival | PLK1 nominally associated with outcome | Exploratory prognostic signal |
| Drug–target | 487 interaction records; 7 genes with associations; 7 multi-target drugs | Translational/hypothesis-generating evidence |

The convergence of these analyses, including the translational drug–target layer, strengthens the biological interpretation while the limitations of sample size, database coverage, and validation scope prevent overstatement of the findings.

---

# 19. Limitations

Several limitations should be considered when interpreting the findings.

1. The discovery cohort contains only 16 samples: 8 Normal and 8 TNBC.
2. The independent validation cohort is also small, containing 3 Normal and 5 TNBC samples.
3. Within-cohort candidate validation does not constitute independent validation.
4. Independent validation reproduced the expression direction of 8/10 genes but not all candidates.
5. The correlation analysis was based on only 8 TNBC samples and is therefore exploratory.
6. The molecular-signature and ROC/AUC analyses were performed in small cohorts, increasing the risk of unstable performance estimates.
7. AUC values of 1.00 should not be interpreted as evidence of clinical diagnostic performance.
8. Network analysis is based on STRING-derived associations and does not establish direct physical interaction or causality.
9. The TCGA-BRCA survival cohort contains 1,073 patients but only 150 observed deaths.
10. Survival analysis used median expression groups and was not adjusted for clinical covariates.
11. PLK1 showed a nominal survival association but did not remain significant after multiple-testing correction.
12. The TCGA-BRCA cohort represents breast cancer broadly and is not equivalent to a TNBC-only cohort unless explicitly restricted.
13. The ten genes should be considered prioritized candidate markers rather than clinically validated biomarkers.
14. Transcriptomic associations do not establish biological causality.
15. DGIdb associations represent database-reported drug–gene relationships and do not establish therapeutic efficacy, drug sensitivity, or clinical suitability.
16. Only 7 of the 10 prioritized genes returned usable drug associations in the analyzed DGIdb query.
17. The number of recorded drug associations reflects database coverage and should not be interpreted as a measure of biological importance or therapeutic potential.
18. Larger, independent, TNBC-specific cohorts and experimental validation are required to establish the reproducibility and biological significance of the identified programs.

---

# 20. Overall Conclusion

The transcriptomic analysis identified extensive molecular differences between Normal and TNBC samples and consistently revealed two major biological programs.

The first was characterized by **increased proliferative and mitotic activity**, supported by differential expression, GO enrichment, KEGG pathways, Hallmark GSEA, coordinated expression of five proliferation-associated candidates, and network-level evidence.

The second was characterized by **reduced lipid-associated and metabolic activity**, supported by differential expression, fatty-acid and metabolic GO terms, PPAR and lipolysis pathways, Hallmark metabolic programs, and coordinated expression of five metabolic-associated candidates.

The final ten-gene panel was:

**CDC20, BUB1, TRIP13, PLK1, AURKB, PNPLA2, PPARG, LIPE, LEP, and CIDEC.**

The independent validation reproduced the expression direction of 8 of these 10 genes, while the two molecular signatures showed significant separation in both the discovery and independent cohorts.

Hallmark GSEA and PPI network analysis provided additional support for the coherence of the two molecular programs. ROC/AUC analysis demonstrated strong separation of the analyzed Normal and TNBC samples, although the small cohorts prevent these values from being interpreted as clinical performance estimates.

In TCGA-BRCA survival analysis, PLK1 showed the strongest nominal association with overall survival, but the association did not remain significant after multiple-testing correction. It should therefore be considered an exploratory survival-associated candidate rather than a validated prognostic biomarker.

DGIdb drug–target analysis added a translational layer by identifying reported drug associations for 7 candidate genes, with the largest drug-association coverage observed for PLK1, PPARG, and AURKB. Seven multi-target compounds connected prioritized genes, including compounds linking PLK1 with PPARG and AURKB with PLK1. These findings are hypothesis-generating and do not establish therapeutic efficacy.

Overall, the findings support a model in which the TNBC transcriptional state examined in this project is characterized by a shift toward **enhanced proliferative activity accompanied by reduced representation of lipid-associated metabolic programs**.

The ten-gene panel provides a biologically coherent starting point for further investigation. Larger TNBC-specific cohorts, independent validation, clinical adjustment, and experimental studies will be necessary to determine whether these candidates have reproducible biomarker or therapeutic relevance.
