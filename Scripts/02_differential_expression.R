# ============================================================
# 02_differential_expression.R
# Differential Expression Analysis using DESeq2
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

# Load DESeq2 for RNA-seq normalization and
# differential expression analysis
library(DESeq2)


# ------------------------------------------------------------
# 2. Prepare metadata for DESeq2
# ------------------------------------------------------------

# Create a copy of the cleaned metadata
# generated during the exploratory analysis
dds_metadata <- analysis_metadata

# Set sample IDs as row names
# DESeq2 uses these to match metadata with count-matrix columns
rownames(dds_metadata) <- dds_metadata$ID_REF

# Convert cancer type into a factor
# Normal is set as the reference group
# Therefore, the analysis represents:
# TNBC compared with Normal
dds_metadata$Cancer_type <- factor(
  dds_metadata$Cancer_type,
  levels = c("Normal", "TNBC")
)

# Convert Patient_ID into a factor
# This allows us to account for patient-to-patient
# variation because the samples are matched within patients
dds_metadata$Patient_ID <- factor(
  dds_metadata$Patient_ID
)

# Remove the ID_REF column because the sample IDs
# are already stored as row names
dds_metadata$ID_REF <- NULL


# ------------------------------------------------------------
# 3. Prepare the count matrix
# ------------------------------------------------------------

# Remove GeneID from the count matrix
# and convert the remaining sample counts into a matrix
count_matrix <- as.matrix(
  counts_subset[, -1]
)

# Set GeneID as the row names
rownames(count_matrix) <- counts_subset$GeneID


# ------------------------------------------------------------
# 4. Create the DESeq2 dataset
# ------------------------------------------------------------

# Create a DESeq2 dataset using the matched-patient design
#
# Patient_ID accounts for differences between patients.
# Cancer_type tests the difference between TNBC and Normal.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = dds_metadata,
  design = ~ Patient_ID + Cancer_type
)

# Display the DESeq2 object
dds


# ------------------------------------------------------------
# 5. Variance-stabilizing transformation
# ------------------------------------------------------------

# Apply variance-stabilizing transformation
# This reduces the dependence of variance on mean expression
# and produces data suitable for visualization
vsd <- vst(
  dds,
  blind = FALSE
)

# Extract the transformed expression matrix
vsd_matrix <- assay(vsd)

# Check the dimensions of the transformed matrix
dim(vsd_matrix)


# ------------------------------------------------------------
# 6. Principal Component Analysis
# ------------------------------------------------------------

# Generate PCA data from the variance-stabilized expression
# using cancer type as the grouping variable
pca_data <- plotPCA(
  vsd,
  intgroup = "Cancer_type",
  returnData = TRUE
)

# Calculate the percentage of variance explained by PC1 and PC2
percent_var <- round(
  100 * attr(pca_data, "percentVar")
)

# Use different symbols for Normal and TNBC
group_symbols <- ifelse(
  pca_data$Cancer_type == "Normal",
  16,
  17
)

# Create the PCA plot
par(mar = c(5, 5, 4, 9))

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

# Add sample labels
text(
  pca_data$PC1,
  pca_data$PC2,
  labels = paste0("S", 1:nrow(pca_data)),
  pos = 3,
  cex = 0.7
)

# Add group legend
legend(
  "topright",
  inset = c(-0.40, 0),
  legend = c("Normal", "TNBC"),
  pch = c(16, 17),
  bty = "n",
  xpd = NA
)

# Reset plotting margins
par(mar = c(5, 4, 4, 2))


# ------------------------------------------------------------
# 7. Run differential expression analysis
# ------------------------------------------------------------

# Run DESeq2
# This performs normalization, dispersion estimation,
# and model fitting using the specified design
dds <- DESeq(dds)

# Extract differential expression results
# The comparison is TNBC relative to Normal
res <- results(
  dds,
  contrast = c(
    "Cancer_type",
    "TNBC",
    "Normal"
  )
)

# Display a summary of the results
summary(res)


# ------------------------------------------------------------
# 8. Convert results into a data frame
# ------------------------------------------------------------

# Convert DESeq2 results into a regular data frame
res_df <- as.data.frame(res)

# Add GeneID as a separate column
res_df$GeneID <- rownames(res_df)

# Display the first few results
head(res_df)

# Sort genes by adjusted p-value
res_df_sorted <- res_df[
  order(res_df$padj),
]

# Display the 20 genes with the smallest adjusted p-values
head(res_df_sorted, 20)


# ------------------------------------------------------------
# 9. Identify differentially expressed genes
# ------------------------------------------------------------

# Define DEGs using:
# adjusted p-value < 0.05
# absolute log2 fold change >= 1
#
# log2FC >= 1  = at least 2-fold higher in TNBC
# log2FC <= -1 = at least 2-fold lower in TNBC

deg <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
]

# Count total significant DEGs
nrow(deg)


# ------------------------------------------------------------
# 10. Separate upregulated and downregulated genes
# ------------------------------------------------------------

# Positive log2FoldChange indicates higher expression in TNBC
upregulated <- deg[
  deg$log2FoldChange >= 1,
]

# Negative log2FoldChange indicates lower expression in TNBC
downregulated <- deg[
  deg$log2FoldChange <= -1,
]

# Count upregulated genes
nrow(upregulated)

# Count downregulated genes
nrow(downregulated)


# ------------------------------------------------------------
# 11. Load gene annotation
# ------------------------------------------------------------

# Load the human gene annotation table provided
# with the GEO dataset
annotation <- read.delim(
  "Data/Human.GRCh38.p13.annot.tsv.gz",
  check.names = FALSE
)

# Keep only the annotation columns required
annotation_small <- annotation[
  ,
  c(
    "GeneID",
    "Symbol",
    "Description",
    "GeneType",
    "EnsemblGeneID"
  )
]

# Add gene annotation to the DEG table
deg_annotated <- merge(
  deg,
  annotation_small,
  by = "GeneID",
  all.x = TRUE
)

# Sort annotated DEGs by adjusted p-value
deg_annotated <- deg_annotated[
  order(deg_annotated$padj),
]

# Display the 20 strongest annotated DEGs
head(deg_annotated, 20)


# ------------------------------------------------------------
# 12. Prepare data for volcano plot
# ------------------------------------------------------------

# Create a copy of the complete DESeq2 results
volcano_data <- res_df

# Remove genes without adjusted p-values
volcano_data <- volcano_data[
  !is.na(volcano_data$padj),
]

# Initially classify all genes as not significant
volcano_data$Significance <- "Not significant"

# Classify significantly upregulated genes
volcano_data$Significance[
  volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange >= 1
] <- "Upregulated"

# Classify significantly downregulated genes
volcano_data$Significance[
  volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange <= -1
] <- "Downregulated"

# Calculate -log10 adjusted p-value
volcano_data$minus_log10_padj <- -log10(
  volcano_data$padj
)

# Display the number of genes in each category
table(volcano_data$Significance)


# ------------------------------------------------------------
# 13. Save differential expression results
# ------------------------------------------------------------

# Create Results folder if it does not already exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save complete DESeq2 results
write.csv(
  res_df,
  "Results/DESeq2_All_Results.csv",
  row.names = FALSE
)

# Save annotated DEG results
write.csv(
  deg_annotated,
  "Results/DEGs_Annotated.csv",
  row.names = FALSE
)

# Save upregulated genes
write.csv(
  upregulated,
  "Results/DEGs_Upregulated.csv",
  row.names = FALSE
)

# Save downregulated genes
write.csv(
  downregulated,
  "Results/DEGs_Downregulated.csv",
  row.names = FALSE
)