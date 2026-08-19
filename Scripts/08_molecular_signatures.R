# ============================================================
# 08_molecular_signatures.R
# Molecular Program / Signature Validation
# ============================================================

# This script evaluates two predefined molecular programs:
#
# 1. Proliferation signature:
#    CDC20, BUB1, TRIP13, PLK1, AURKB
#
# 2. Metabolic signature:
#    PNPLA2, PPARG, LIPE, LEP, CIDEC
#
# The signatures are evaluated in:
# - The original discovery cohort
# - The independent GSE52194 validation cohort
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(DESeq2)


# ------------------------------------------------------------
# 2. Create output directories
# ------------------------------------------------------------

if (!dir.exists("Results/Molecular_Signatures")) {
  dir.create(
    "Results/Molecular_Signatures",
    recursive = TRUE
  )
}

if (!dir.exists("Figures/Molecular_Signatures")) {
  dir.create(
    "Figures/Molecular_Signatures",
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Load objects from previous analyses
# ------------------------------------------------------------

# Script 01 creates the original discovery-cohort
# expression objects used for signature construction.
source("Scripts/01_explore_data.R")

# Script 07 performs the independent GSE52194 validation
# and creates the independent-cohort expression objects.
source("Scripts/07_independent_validation.R")


# ------------------------------------------------------------
# 4. Define the two molecular signatures
# ------------------------------------------------------------

# Proliferation-associated genes identified in the
# original TNBC analysis.
proliferation_genes <- c(
  "CDC20",
  "BUB1",
  "TRIP13",
  "PLK1",
  "AURKB"
)

# Metabolic/lipid-associated genes identified in the
# original TNBC analysis.
metabolic_genes <- c(
  "PNPLA2",
  "PPARG",
  "LIPE",
  "LEP",
  "CIDEC"
)


# ------------------------------------------------------------
# 5. Map candidate genes to the original expression matrix
# ------------------------------------------------------------

# The original expression matrix uses numeric GeneIDs,
# while the signatures are defined using gene symbols.
candidate_annotation <- annotation[
  annotation$Symbol %in% c(
    proliferation_genes,
    metabolic_genes
  ),
  c(
    "GeneID",
    "Symbol",
    "EnsemblGeneID"
  ),
  drop = FALSE
]

# Verify that all ten candidate genes were mapped.
if (nrow(candidate_annotation) != 10) {
  stop(
    "Expected 10 candidate genes, but found ",
    nrow(candidate_annotation),
    "."
  )
}


# ------------------------------------------------------------
# 6. Extract candidate-gene expression
# ------------------------------------------------------------

candidate_gene_ids_original <-
  candidate_annotation$GeneID

candidate_vst_expression <- vsd_matrix[
  as.character(candidate_gene_ids_original),
  ,
  drop = FALSE
]

# Replace numeric GeneIDs with gene symbols.
rownames(candidate_vst_expression) <-
  candidate_annotation$Symbol



# ------------------------------------------------------------
# 7. Calculate discovery-cohort signature scores
# ------------------------------------------------------------

# Standardize each gene across the discovery samples.
candidate_z_expression <- t(
  scale(
    t(candidate_vst_expression)
  )
)

# Calculate the mean standardized expression of the
# five proliferation-associated genes.
proliferation_score <- colMeans(
  candidate_z_expression[
    proliferation_genes,
    ,
    drop = FALSE
  ]
)

# Calculate the mean standardized expression of the
# five metabolic-associated genes.
metabolic_score <- colMeans(
  candidate_z_expression[
    metabolic_genes,
    ,
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 8. Create discovery-cohort signature table
# ------------------------------------------------------------

signature_scores <- data.frame(
  Sample = names(proliferation_score),
  Group = colData(dds)$Cancer_type[
    match(
      names(proliferation_score),
      rownames(colData(dds))
    )
  ],
  Proliferation_Score = as.numeric(
    proliferation_score
  ),
  Metabolic_Score = as.numeric(
    metabolic_score
  ),
  stringsAsFactors = FALSE
)

signature_scores$Group <- factor(
  signature_scores$Group,
  levels = c("Normal", "TNBC")
)


# ------------------------------------------------------------
# 9. Statistical comparison of discovery-cohort signatures
# ------------------------------------------------------------

# Test whether the proliferation signature differs
# between Normal and TNBC samples.
proliferation_test <- wilcox.test(
  Proliferation_Score ~ Group,
  data = signature_scores,
  exact = FALSE
)

# Test whether the metabolic signature differs
# between Normal and TNBC samples.
metabolic_test <- wilcox.test(
  Metabolic_Score ~ Group,
  data = signature_scores,
  exact = FALSE
)


# ------------------------------------------------------------
# 10. Summarize discovery-cohort signature differences
# ------------------------------------------------------------

signature_summary <- aggregate(
  cbind(
    Proliferation_Score,
    Metabolic_Score
  ) ~ Group,
  data = signature_scores,
  FUN = mean
)

signature_medians <- aggregate(
  cbind(
    Proliferation_Score,
    Metabolic_Score
  ) ~ Group,
  data = signature_scores,
  FUN = median
)


# ------------------------------------------------------------
# 11. Extract signatures from independent GSE52194 cohort
# ------------------------------------------------------------

validation_signature_expression <-
  candidate_normalized_expression[
    c(
      proliferation_genes,
      metabolic_genes
    ),
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 12. Calculate independent-cohort signature scores
# ------------------------------------------------------------

# Standardize each candidate gene across the independent
# validation samples.
validation_z_expression <- t(
  scale(
    t(validation_signature_expression)
  )
)

# Calculate the proliferation signature score.
validation_proliferation_score <- colMeans(
  validation_z_expression[
    proliferation_genes,
    ,
    drop = FALSE
  ]
)

# Calculate the metabolic signature score.
validation_metabolic_score <- colMeans(
  validation_z_expression[
    metabolic_genes,
    ,
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 13. Create independent-cohort signature table
# ------------------------------------------------------------

validation_signature_scores <- data.frame(
  Sample = names(validation_proliferation_score),
  Group = validation_coldata$Group[
    match(
      names(validation_proliferation_score),
      rownames(validation_coldata)
    )
  ],
  Proliferation_Score = as.numeric(
    validation_proliferation_score
  ),
  Metabolic_Score = as.numeric(
    validation_metabolic_score
  ),
  stringsAsFactors = FALSE
)

validation_signature_scores$Group <- factor(
  validation_signature_scores$Group,
  levels = c("Normal", "TNBC")
)


# ------------------------------------------------------------
# 14. Statistical comparison in GSE52194
# ------------------------------------------------------------

validation_proliferation_test <- wilcox.test(
  Proliferation_Score ~ Group,
  data = validation_signature_scores,
  exact = FALSE
)

validation_metabolic_test <- wilcox.test(
  Metabolic_Score ~ Group,
  data = validation_signature_scores,
  exact = FALSE
)


# ------------------------------------------------------------
# 15. Summarize independent-cohort signature differences
# ------------------------------------------------------------

validation_signature_summary <- aggregate(
  cbind(
    Proliferation_Score,
    Metabolic_Score
  ) ~ Group,
  data = validation_signature_scores,
  FUN = mean
)

validation_signature_medians <- aggregate(
  cbind(
    Proliferation_Score,
    Metabolic_Score
  ) ~ Group,
  data = validation_signature_scores,
  FUN = median
)


# ------------------------------------------------------------
# 16. Save signature scores
# ------------------------------------------------------------

write.csv(
  signature_scores,
  "Results/Molecular_Signatures/Discovery_Cohort_Signature_Scores.csv",
  row.names = FALSE
)

write.csv(
  validation_signature_scores,
  "Results/Molecular_Signatures/GSE52194_Signature_Scores.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Create discovery-cohort signature figure
# ------------------------------------------------------------

png(
  "Figures/Molecular_Signatures/Discovery_Cohort_Molecular_Signatures.png",
  width = 2400,
  height = 1800,
  res = 250
)

par(
  mfrow = c(1, 2),
  mar = c(5, 5, 4, 1)
)

boxplot(
  Proliferation_Score ~ Group,
  data = signature_scores,
  main = "Proliferation Signature",
  xlab = "",
  ylab = "Signature Score"
)

stripchart(
  Proliferation_Score ~ Group,
  data = signature_scores,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  add = TRUE
)

boxplot(
  Metabolic_Score ~ Group,
  data = signature_scores,
  main = "Metabolic Signature",
  xlab = "",
  ylab = "Signature Score"
)

stripchart(
  Metabolic_Score ~ Group,
  data = signature_scores,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  add = TRUE
)

dev.off()

par(
  mfrow = c(1, 1)
)


# ------------------------------------------------------------
# 18. Create independent-cohort signature figure
# ------------------------------------------------------------

png(
  "Figures/Molecular_Signatures/GSE52194_Molecular_Signatures.png",
  width = 2400,
  height = 1800,
  res = 250
)

par(
  mfrow = c(1, 2),
  mar = c(5, 5, 4, 1)
)

boxplot(
  Proliferation_Score ~ Group,
  data = validation_signature_scores,
  main = "Proliferation Signature",
  xlab = "",
  ylab = "Signature Score"
)

stripchart(
  Proliferation_Score ~ Group,
  data = validation_signature_scores,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  add = TRUE
)

boxplot(
  Metabolic_Score ~ Group,
  data = validation_signature_scores,
  main = "Metabolic Signature",
  xlab = "",
  ylab = "Signature Score"
)

stripchart(
  Metabolic_Score ~ Group,
  data = validation_signature_scores,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  add = TRUE
)

dev.off()

par(
  mfrow = c(1, 1)
)


# ------------------------------------------------------------
# 19. Print summary
# ------------------------------------------------------------

cat(
  "\nMolecular signature analysis completed.\n"
)

cat(
  "Discovery cohort: 8 Normal, 8 TNBC\n"
)

cat(
  "Discovery proliferation p-value:",
  proliferation_test$p.value,
  "\n"
)

cat(
  "Discovery metabolic p-value:",
  metabolic_test$p.value,
  "\n"
)

cat(
  "Independent cohort: 3 Normal, 5 TNBC\n"
)

cat(
  "Independent proliferation p-value:",
  validation_proliferation_test$p.value,
  "\n"
)

cat(
  "Independent metabolic p-value:",
  validation_metabolic_test$p.value,
  "\n"
)