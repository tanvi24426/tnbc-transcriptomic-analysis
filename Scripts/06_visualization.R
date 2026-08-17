# ============================================================
# 06_visualization.R
# Visualization of Differential Expression and
# Final Candidate Gene Results
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(DESeq2)
library(pheatmap)


# ------------------------------------------------------------
# 2. Create the Figures folder
# ------------------------------------------------------------

# Create a dedicated folder for all final figures
if (!dir.exists("Figures")) {
  dir.create("Figures")
}


# ------------------------------------------------------------
# 3. PCA plot
# ------------------------------------------------------------

# PCA was calculated using the variance-stabilized expression
# data generated during differential expression analysis.

pca_data <- plotPCA(
  vsd,
  intgroup = "Cancer_type",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar")
)

# Define plotting symbols for each sample group
group_symbols <- ifelse(
  pca_data$Cancer_type == "Normal",
  16,
  17
)

# Save PCA plot as a PNG file
png(
  "Figures/PCA_Normal_vs_TNBC.png",
  width = 1800,
  height = 1600,
  res = 250
)

par(
  mar = c(5, 5, 4, 9)
)

plot(
  pca_data$PC1,
  pca_data$PC2,
  pch = group_symbols,
  xlab = paste0(
    "PC1: ",
    percent_var[1],
    "% variance"
  ),
  ylab = paste0(
    "PC2: ",
    percent_var[2],
    "% variance"
  ),
  main = "PCA of Normal and TNBC Samples"
)

text(
  pca_data$PC1,
  pca_data$PC2,
  labels = paste0(
    "S",
    1:nrow(pca_data)
  ),
  pos = 3,
  cex = 0.7
)

legend(
  "topright",
  inset = c(-0.40, 0),
  legend = c(
    "Normal",
    "TNBC"
  ),
  pch = c(16, 17),
  bty = "n",
  xpd = NA
)

dev.off()

par(
  mar = c(5, 4, 4, 2)
)


# ------------------------------------------------------------
# 4. Volcano plot
# ------------------------------------------------------------

# Prepare the volcano plot data
volcano_data <- res_df[
  !is.na(res_df$padj),
]

# Assign significance categories
volcano_data$Significance <- "Not significant"

volcano_data$Significance[
  volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange >= 1
] <- "Upregulated"

volcano_data$Significance[
  volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange <= -1
] <- "Downregulated"

# Calculate negative log10 adjusted p-value
volcano_data$minus_log10_padj <- -log10(
  volcano_data$padj
)

# Assign plotting symbols
# 17 = Upregulated
# 15 = Downregulated
# 16 = Not significant
volcano_symbols <- ifelse(
  volcano_data$Significance == "Upregulated",
  17,
  ifelse(
    volcano_data$Significance == "Downregulated",
    15,
    16
  )
)

# Save volcano plot
png(
  "Figures/Volcano_Plot_TNBC_vs_Normal.png",
  width = 2200,
  height = 2000,
  res = 250
)

plot(
  volcano_data$log2FoldChange,
  volcano_data$minus_log10_padj,
  pch = volcano_symbols,
  cex = 0.55,
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Differential Expression: TNBC vs Normal"
)

# Add significance threshold lines
abline(
  v = c(-1, 1),
  lty = 2
)

abline(
  h = -log10(0.05),
  lty = 2
)

# Identify the final prioritized genes
label_indices <- which(
  rownames(volcano_data) %in% final_gene_ids
)

# Add labels for the 10 prioritized genes
if (length(label_indices) > 0) {
  
  label_symbols <- final_gene_symbols[
    match(
      rownames(volcano_data)[label_indices],
      final_gene_ids
    )
  ]
  
  text(
    volcano_data$log2FoldChange[label_indices],
    volcano_data$minus_log10_padj[label_indices],
    labels = label_symbols,
    pos = 3,
    cex = 0.7
  )
}

# Add legend
legend(
  "topright",
  legend = c(
    "Upregulated",
    "Downregulated",
    "Not significant"
  ),
  pch = c(
    17,
    15,
    16
  ),
  bty = "n"
)

dev.off()


# ------------------------------------------------------------
# 5. Final 10-gene expression heatmap
# ------------------------------------------------------------

# Extract the final 10 genes from the VST matrix
candidate_expression <- vsd_matrix[
  final_gene_ids,
  ,
  drop = FALSE
]

# Replace GeneIDs with gene symbols
rownames(candidate_expression) <- final_gene_symbols

# Create sample annotation
sample_annotation <- data.frame(
  Cancer_type = dds_metadata[
    colnames(candidate_expression),
    "Cancer_type"
  ]
)

rownames(sample_annotation) <- colnames(
  candidate_expression
)

# Save the final 10-gene heatmap
png(
  "Figures/Final_10_Gene_Heatmap.png",
  width = 2200,
  height = 2000,
  res = 250
)

pheatmap(
  candidate_expression,
  scale = "row",
  annotation_col = sample_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Expression of Prioritized TNBC Candidate Genes",
  fontsize_row = 10,
  fontsize_col = 8
)

dev.off()


# ------------------------------------------------------------
# 6. Final 10-gene boxplots
# ------------------------------------------------------------

# Convert the expression matrix into long-format data
boxplot_data <- candidate_expression_df

# Create one figure containing all candidate-gene comparisons
png(
  "Figures/Final_10_Gene_Boxplots.png",
  width = 2400,
  height = 2200,
  res = 250
)

par(
  mfrow = c(5, 2),
  mar = c(5, 4, 3, 1)
)

for (gene in final_gene_symbols) {
  
  # Extract expression values for the current gene
  expression_values <- boxplot_data[
    ,
    gene
  ]
  
  # Create boxplot by cancer type
  boxplot(
    expression_values ~ boxplot_data$Cancer_type,
    main = gene,
    xlab = "Cancer Type",
    ylab = "VST Expression"
  )
  
  # Add individual sample points
  stripchart(
    expression_values ~ boxplot_data$Cancer_type,
    vertical = TRUE,
    method = "jitter",
    pch = 16,
    add = TRUE
  )
}

dev.off()

# Reset plotting layout
par(
  mfrow = c(1, 1)
)


# ------------------------------------------------------------
# 7. TNBC gene-gene correlation heatmap
# ------------------------------------------------------------

# Extract only TNBC samples
tnbc_expression <- candidate_expression[
  ,
  tnbc_samples,
  drop = FALSE
]

# Calculate Pearson correlations
candidate_correlations <- cor(
  t(tnbc_expression),
  method = "pearson"
)

# Save correlation heatmap
png(
  "Figures/TNBC_10_Gene_Correlation_Heatmap.png",
  width = 2200,
  height = 2000,
  res = 250
)

pheatmap(
  candidate_correlations,
  color = colorRampPalette(
    c(
      "blue",
      "white",
      "red"
    )
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
  main = "Gene-Gene Correlation in TNBC",
  fontsize_row = 9,
  fontsize_col = 9
)

dev.off()


# ------------------------------------------------------------
# 8. Confirm generated figures
# ------------------------------------------------------------

# List all PNG figures generated by this script
list.files(
  "Figures",
  pattern = "\\.png$"
)