# ============================================================
# 07_independent_validation.R
# Independent Validation of the 10-Gene TNBC Candidate Panel
# Dataset: GSE52194
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(DESeq2)


# ------------------------------------------------------------
# 2. Create output directories
# ------------------------------------------------------------

if (!dir.exists("Results/Independent_Validation")) {
  dir.create(
    "Results/Independent_Validation",
    recursive = TRUE
  )
}

if (!dir.exists("Figures/Independent_Validation")) {
  dir.create(
    "Figures/Independent_Validation",
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Load GSE52194 raw-count data
# ------------------------------------------------------------

validation_file <-
  "Data/Independent_Validation/E-GEOD-52194-raw-counts.tsv"

validation_counts <- read.delim(
  validation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Load experiment design
# ------------------------------------------------------------

validation_design <- read.delim(
  "https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-52194/resources/ExperimentDesignFile.RnaSeq/experiment-design",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Prepare sample metadata
# ------------------------------------------------------------

validation_metadata <- validation_design[
  validation_design$Run %in% colnames(validation_counts),
  c(
    "Run",
    "Sample Characteristic[clinical information]",
    "Sample Characteristic[individual]",
    "Sample Characteristic[organism part]",
    "Analysed"
  )
]

colnames(validation_metadata) <- c(
  "Run",
  "Group",
  "Individual",
  "Organism_part",
  "Analysed"
)


# ------------------------------------------------------------
# 6. Select analyzed TNBC and Normal samples
# ------------------------------------------------------------

validation_metadata <- validation_metadata[
  validation_metadata$Group %in% c(
    "triple-negative breast cancer",
    "normal"
  ) &
    validation_metadata$Analysed == "Yes",
]

validation_sample_ids <- validation_metadata$Run


# ------------------------------------------------------------
# 7. Create validation count matrix
# ------------------------------------------------------------

validation_matrix <- validation_counts[
  ,
  c("Gene ID", "Gene Name", validation_sample_ids),
  drop = FALSE
]

validation_count_matrix <- as.matrix(
  validation_matrix[
    ,
    validation_sample_ids,
    drop = FALSE
  ]
)

# Use Ensembl Gene IDs as row names.
rownames(validation_count_matrix) <-
  validation_matrix$`Gene ID`


# ------------------------------------------------------------
# 8. Create validation sample metadata
# ------------------------------------------------------------

validation_coldata <- data.frame(
  Group = factor(
    ifelse(
      validation_metadata$Group ==
        "triple-negative breast cancer",
      "TNBC",
      "Normal"
    ),
    levels = c("Normal", "TNBC")
  )
)

rownames(validation_coldata) <-
  validation_metadata$Run


# ------------------------------------------------------------
# 9. Normalize validation counts
# ------------------------------------------------------------

validation_dds <- DESeqDataSetFromMatrix(
  countData = round(validation_count_matrix),
  colData = validation_coldata,
  design = ~ Group
)

validation_dds <- estimateSizeFactors(
  validation_dds
)

validation_normalized_counts <- counts(
  validation_dds,
  normalized = TRUE
)


# ------------------------------------------------------------
# 10. Define the frozen 10-gene candidate panel
# ------------------------------------------------------------

final_validation_genes <- c(
  "CDC20",
  "BUB1",
  "TRIP13",
  "PLK1",
  "AURKB",
  "PNPLA2",
  "PPARG",
  "LIPE",
  "LEP",
  "CIDEC"
)

candidate_gene_ids <- c(
  CDC20  = "ENSG00000117399",
  BUB1   = "ENSG00000169679",
  TRIP13 = "ENSG00000071539",
  PLK1   = "ENSG00000166851",
  AURKB  = "ENSG00000178999",
  PNPLA2 = "ENSG00000177666",
  PPARG  = "ENSG00000132170",
  LIPE   = "ENSG00000079435",
  LEP    = "ENSG00000174697",
  CIDEC  = "ENSG00000187288"
)


# ------------------------------------------------------------
# 11. Extract normalized expression
# ------------------------------------------------------------

candidate_normalized_expression <-
  validation_normalized_counts[
    candidate_gene_ids,
    ,
    drop = FALSE
  ]

rownames(candidate_normalized_expression) <-
  names(candidate_gene_ids)


# ------------------------------------------------------------
# 12. Calculate validation statistics
# ------------------------------------------------------------

validation_results <- data.frame(
  Gene = final_validation_genes,
  Normal_Mean = NA_real_,
  TNBC_Mean = NA_real_,
  Log2_Fold_Change = NA_real_,
  P_Value = NA_real_,
  stringsAsFactors = FALSE
)

for (gene in final_validation_genes) {
  
  expression_values <- as.numeric(
    candidate_normalized_expression[gene, ]
  )
  
  normal_values <- expression_values[
    validation_coldata$Group == "Normal"
  ]
  
  tnbc_values <- expression_values[
    validation_coldata$Group == "TNBC"
  ]
  
  validation_results[
    validation_results$Gene == gene,
    "Normal_Mean"
  ] <- mean(normal_values)
  
  validation_results[
    validation_results$Gene == gene,
    "TNBC_Mean"
  ] <- mean(tnbc_values)
  
  validation_results[
    validation_results$Gene == gene,
    "Log2_Fold_Change"
  ] <- log2(
    (mean(tnbc_values) + 1) /
      (mean(normal_values) + 1)
  )
  
  validation_results[
    validation_results$Gene == gene,
    "P_Value"
  ] <- wilcox.test(
    tnbc_values,
    normal_values,
    exact = FALSE
  )$p.value
}


# ------------------------------------------------------------
# 13. Multiple-testing correction
# ------------------------------------------------------------

validation_results$Adjusted_P_Value <- p.adjust(
  validation_results$P_Value,
  method = "BH"
)


# ------------------------------------------------------------
# 14. Determine validation direction
# ------------------------------------------------------------

validation_results$Validation_Direction <- ifelse(
  validation_results$Log2_Fold_Change > 0,
  "Up in TNBC",
  "Down in TNBC"
)


# ------------------------------------------------------------
# 15. Define original direction
# ------------------------------------------------------------

validation_results$Original_Direction <- ifelse(
  validation_results$Gene %in% c(
    "CDC20",
    "BUB1",
    "TRIP13",
    "PLK1",
    "AURKB"
  ),
  "Up in TNBC",
  "Down in TNBC"
)


# ------------------------------------------------------------
# 16. Determine directional replication
# ------------------------------------------------------------

validation_results$Direction_Replicated <-
  validation_results$Validation_Direction ==
  validation_results$Original_Direction


# ------------------------------------------------------------
# 17. Save validation results
# ------------------------------------------------------------

write.csv(
  validation_results,
  "Results/Independent_Validation/GSE52194_10_Gene_Validation.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Save replication summary
# ------------------------------------------------------------

replication_summary_df <- data.frame(
  Result = c(
    "Direction reproduced",
    "Direction not reproduced"
  ),
  Number_of_Genes = c(
    sum(validation_results$Direction_Replicated),
    sum(!validation_results$Direction_Replicated)
  )
)

write.csv(
  replication_summary_df,
  "Results/Independent_Validation/GSE52194_Replication_Summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Create validation plotting data
# ------------------------------------------------------------

validation_plot_data <- data.frame(
  Gene = rep(
    rownames(candidate_normalized_expression),
    each = ncol(candidate_normalized_expression)
  ),
  Sample = rep(
    colnames(candidate_normalized_expression),
    times = nrow(candidate_normalized_expression)
  ),
  Expression = as.vector(
    t(candidate_normalized_expression)
  )
)

validation_plot_data$Group <- validation_coldata[
  validation_plot_data$Sample,
  "Group"
]


# ------------------------------------------------------------
# 20. Generate validation figure
# ------------------------------------------------------------

png(
  "Figures/Independent_Validation/GSE52194_10_Gene_Expression_Validation.png",
  width = 2600,
  height = 2200,
  res = 250
)

par(
  mfrow = c(5, 2),
  mar = c(5, 4, 3, 1)
)

for (gene in final_validation_genes) {
  
  expression_values <- as.numeric(
    candidate_normalized_expression[gene, ]
  )
  
  normal_values <- expression_values[
    validation_coldata$Group == "Normal"
  ]
  
  tnbc_values <- expression_values[
    validation_coldata$Group == "TNBC"
  ]
  
  boxplot(
    normal_values,
    tnbc_values,
    names = c("Normal", "TNBC"),
    main = gene,
    ylab = "Normalized Expression",
    xlab = "Group"
  )
  
  stripchart(
    list(
      normal_values,
      tnbc_values
    ),
    vertical = TRUE,
    method = "jitter",
    pch = 16,
    add = TRUE
  )
}

dev.off()

par(
  mfrow = c(1, 1)
)


# ------------------------------------------------------------
# 21. Print validation summary
# ------------------------------------------------------------

cat(
  "\nIndependent validation completed.\n"
)

cat(
  "Normal samples:",
  sum(validation_coldata$Group == "Normal"),
  "\n"
)

cat(
  "TNBC samples:",
  sum(validation_coldata$Group == "TNBC"),
  "\n"
)

cat(
  "Genes reproducing original direction:",
  sum(validation_results$Direction_Replicated),
  "of",
  length(final_validation_genes),
  "\n"
)