# ============================================================
# 05_candidate_validation.R
# Validation of the Final 10 Candidate Genes
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(pheatmap)


# ------------------------------------------------------------
# 2. Define the final 10 candidate genes
# ------------------------------------------------------------

final_gene_symbols <- c(
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


# ------------------------------------------------------------
# 3. Extract expression values for final candidates
# ------------------------------------------------------------

# Match the final candidate gene symbols to their GeneIDs
final_gene_ids <- final_candidates$GeneID[
  match(
    final_gene_symbols,
    final_candidates$Symbol
  )
]

# Extract variance-stabilized expression values
# for the final 10 genes
candidate_expression <- vsd_matrix[
  final_gene_ids,
  ,
  drop = FALSE
]

# Replace GeneIDs with gene symbols for easier interpretation
rownames(candidate_expression) <- final_gene_symbols

# Display the expression matrix
candidate_expression


# ------------------------------------------------------------
# 4. Create a sample-level expression table
# ------------------------------------------------------------

# Transpose the expression matrix so that:
# rows = samples
# columns = genes
candidate_expression_df <- as.data.frame(
  t(candidate_expression)
)

# Add sample IDs
candidate_expression_df$Sample <- rownames(
  candidate_expression_df
)

# Add cancer type information
candidate_expression_df$Cancer_type <- dds_metadata[
  rownames(candidate_expression_df),
  "Cancer_type"
]

# Display the first few samples
head(candidate_expression_df)


# ------------------------------------------------------------
# 5. Identify TNBC samples
# ------------------------------------------------------------

# Extract sample IDs belonging to the TNBC group
tnbc_samples <- rownames(
  candidate_expression_df[
    candidate_expression_df$Cancer_type == "TNBC",
  ]
)

# Display TNBC sample IDs
tnbc_samples

# Count TNBC samples
length(tnbc_samples)


# ------------------------------------------------------------
# 6. Calculate gene–gene correlations within TNBC
# ------------------------------------------------------------

# Extract expression values using only TNBC samples
tnbc_expression <- candidate_expression[
  ,
  tnbc_samples,
  drop = FALSE
]

# Calculate pairwise Pearson correlations
# between the 10 candidate genes across TNBC samples
candidate_correlations <- cor(
  t(tnbc_expression),
  method = "pearson"
)

# Display the correlation matrix
round(
  candidate_correlations,
  2
)


# ------------------------------------------------------------
# 7. Generate TNBC correlation heatmap
# ------------------------------------------------------------

# Create a heatmap showing gene–gene correlations
# among the 10 prioritized candidates in TNBC samples

pheatmap(
  candidate_correlations,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  breaks = seq(
    -1,
    1,
    length.out = 101
  ),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Gene–Gene Correlation in TNBC",
  fontsize_row = 9,
  fontsize_col = 9
)


# ------------------------------------------------------------
# 8. Statistical validation between Normal and TNBC
# ------------------------------------------------------------

# Create an empty list to store the results
validation_results <- list()

# Perform a Wilcoxon rank-sum test for each candidate gene
for (gene in final_gene_symbols) {
  
  # Extract expression values for the current gene
  gene_expression <- candidate_expression_df[
    ,
    gene
  ]
  
  # Separate Normal and TNBC expression values
  normal_values <- gene_expression[
    candidate_expression_df$Cancer_type == "Normal"
  ]
  
  tnbc_values <- gene_expression[
    candidate_expression_df$Cancer_type == "TNBC"
  ]
  
  # Calculate mean expression in each group
  normal_mean <- mean(
    normal_values
  )
  
  tnbc_mean <- mean(
    tnbc_values
  )
  
  # Calculate the difference between groups
  mean_difference <- tnbc_mean - normal_mean
  
  # Perform Wilcoxon rank-sum test
  wilcox_result <- wilcox.test(
    normal_values,
    tnbc_values,
    exact = FALSE
  )
  
  # Store results for the current gene
  validation_results[[gene]] <- data.frame(
    Symbol = gene,
    Normal_Mean = normal_mean,
    TNBC_Mean = tnbc_mean,
    Mean_Difference = mean_difference,
    Wilcoxon_pvalue = wilcox_result$p.value
  )
}


# ------------------------------------------------------------
# 9. Combine statistical validation results
# ------------------------------------------------------------

# Combine individual gene results
validation_table <- do.call(
  rbind,
  validation_results
)

# Adjust p-values using the Benjamini-Hochberg method
validation_table$Adjusted_pvalue <- p.adjust(
  validation_table$Wilcoxon_pvalue,
  method = "BH"
)

# Display the validation results
validation_table


# ------------------------------------------------------------
# 10. Save validation results
# ------------------------------------------------------------

# Create Results folder if necessary
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save the final statistical validation table
write.csv(
  validation_table,
  "Results/Final_10_Gene_Statistical_Validation.csv",
  row.names = FALSE
)

# Save the TNBC correlation matrix
write.csv(
  candidate_correlations,
  "Results/TNBC_10_Gene_Correlation_Matrix.csv",
  row.names = TRUE
)
