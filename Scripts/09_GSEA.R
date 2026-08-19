# ============================================================
# 09_GSEA.R
# Hallmark Gene Set Enrichment Analysis
# ============================================================

# ============================================================
# 1. Load required packages
# ============================================================

library(fgsea)


# ============================================================
# 2. Prepare ranked gene list
# ============================================================

# Use the DESeq2 test statistic as the ranking metric.
# Positive values indicate genes shifted toward TNBC,
# while negative values indicate genes shifted toward Normal.

gsea_ranking <- merge(
  res_df[, c("GeneID", "stat")],
  annotation[, c("GeneID", "Symbol", "EnsemblGeneID")],
  by = "GeneID"
)

# Keep genes with valid statistics and Ensembl identifiers.
gsea_ranking <- gsea_ranking[
  !is.na(gsea_ranking$stat) &
    !is.na(gsea_ranking$EnsemblGeneID) &
    gsea_ranking$EnsemblGeneID != "",
  ,
  drop = FALSE
]

# Remove duplicate Ensembl identifiers.
gsea_ranking <- gsea_ranking[
  !duplicated(gsea_ranking$EnsemblGeneID),
  ,
  drop = FALSE
]

# Sort by DESeq2 test statistic.
gsea_ranking <- gsea_ranking[
  order(
    gsea_ranking$stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


# ============================================================
# 3. Create gene-symbol ranked list
# ============================================================

# Hallmark GMT files use gene symbols.

gsea_symbol_ranking <- gsea_ranking[
  !is.na(gsea_ranking$Symbol) &
    gsea_ranking$Symbol != "",
  ,
  drop = FALSE
]

# Remove duplicate gene symbols.
gsea_symbol_ranking <- gsea_symbol_ranking[
  !duplicated(gsea_symbol_ranking$Symbol),
  ,
  drop = FALSE
]

# Re-sort after removing duplicates.
gsea_symbol_ranking <- gsea_symbol_ranking[
  order(
    gsea_symbol_ranking$stat,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

# Create the ranked vector required by fgsea.
gsea_symbol_list <- gsea_symbol_ranking$stat

names(gsea_symbol_list) <- gsea_symbol_ranking$Symbol


# ============================================================
# 4. Obtain Hallmark gene sets
# ============================================================

# Official MSigDB 2026.1.Hs Hallmark gene sets.
hallmark_url <- paste0(
  "https://data.broadinstitute.org/",
  "gsea-msigdb/msigdb/release/2026.1.Hs/",
  "h.all.v2026.1.Hs.symbols.gmt"
)

# Create local directory for the gene-set file.
if (!dir.exists("Data/GSEA")) {
  dir.create(
    "Data/GSEA",
    recursive = TRUE
  )
}

hallmark_file <- paste0(
  "Data/GSEA/",
  "h.all.v2026.1.Hs.symbols.gmt"
)

# Download only if the file is not already present.
if (!file.exists(hallmark_file)) {
  
  download.file(
    hallmark_url,
    hallmark_file,
    mode = "wb"
  )
  
}

# Confirm that the file exists.
if (!file.exists(hallmark_file)) {
  stop(
    "Hallmark GMT file could not be found."
  )
}

# Read Hallmark pathways.
hallmark_pathways <- gmtPathways(
  hallmark_file
)


# ============================================================
# 5. Run Hallmark GSEA
# ============================================================

# GSEA uses the complete ranked gene list.
#
# Ties:
# Approximately 4.26% of the DESeq2 statistics were tied
# Set a seed for reproducible fgsea results.
set.seed(12345)

# Run fgsea using a single process to improve reproducibility
# and avoid parallel serialization warnings.
gsea_results <- fgsea(
  pathways = hallmark_pathways,
  stats = gsea_symbol_list,
  minSize = 15,
  maxSize = 500,
  eps = 0,
  nproc = 1
)

# Convert to data frame.
gsea_results <- as.data.frame(
  gsea_results
)

# Sort by adjusted p-value.
gsea_results <- gsea_results[
  order(
    gsea_results$padj
  ),
  ,
  drop = FALSE
]


# ============================================================
# 6. Identify significant pathways
# ============================================================

gsea_significant <- gsea_results[
  gsea_results$padj < 0.05,
  ,
  drop = FALSE
]

# Positive NES = TNBC-enriched.
gsea_tnbc <- gsea_significant[
  gsea_significant$NES > 0,
  ,
  drop = FALSE
]

# Negative NES = Normal-enriched.
gsea_normal <- gsea_significant[
  gsea_significant$NES < 0,
  ,
  drop = FALSE
]


# ============================================================
# 7. Prepare results for CSV export
# ============================================================

# fgsea contains a list-column called leadingEdge.
# Convert it to a semicolon-separated character string
# so the results can be exported as CSV.

gsea_results_export <- gsea_results

gsea_results_export$leadingEdge <- vapply(
  gsea_results_export$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

gsea_significant_export <- gsea_significant

gsea_significant_export$leadingEdge <- vapply(
  gsea_significant_export$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

gsea_tnbc_export <- gsea_tnbc

gsea_tnbc_export$leadingEdge <- vapply(
  gsea_tnbc_export$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

gsea_normal_export <- gsea_normal

gsea_normal_export$leadingEdge <- vapply(
  gsea_normal_export$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)


# ============================================================
# 8. Save GSEA results
# ============================================================

if (!dir.exists("Results/GSEA")) {
  dir.create(
    "Results/GSEA",
    recursive = TRUE
  )
}

write.csv(
  gsea_results_export,
  "Results/GSEA/Hallmark_GSEA_All_Results.csv",
  row.names = FALSE
)

write.csv(
  gsea_significant_export,
  "Results/GSEA/Hallmark_GSEA_Significant.csv",
  row.names = FALSE
)

write.csv(
  gsea_tnbc_export,
  "Results/GSEA/Hallmark_GSEA_TNBC_Enriched.csv",
  row.names = FALSE
)

write.csv(
  gsea_normal_export,
  "Results/GSEA/Hallmark_GSEA_Normal_Enriched.csv",
  row.names = FALSE
)


# ============================================================
# 9. Select pathways for visualization
# ============================================================

top_tnbc <- head(
  gsea_tnbc[
    order(gsea_tnbc$padj),
    ,
    drop = FALSE
  ],
  8
)

top_normal <- head(
  gsea_normal[
    order(gsea_normal$padj),
    ,
    drop = FALSE
  ],
  8
)

gsea_plot_data <- rbind(
  top_tnbc,
  top_normal
)

gsea_plot_data <- gsea_plot_data[
  order(gsea_plot_data$NES),
  ,
  drop = FALSE
]


# ============================================================
# 10. Prepare pathway labels
# ============================================================

gsea_plot_data$Pathway <- gsub(
  "^HALLMARK_",
  "",
  gsea_plot_data$pathway
)

gsea_plot_data$Pathway <- gsub(
  "_",
  " ",
  gsea_plot_data$Pathway
)

gsea_plot_data$Pathway <- factor(
  gsea_plot_data$Pathway,
  levels = gsea_plot_data$Pathway
)


# ============================================================
# 11. Create GSEA figure
# ============================================================

if (!dir.exists("Figures/GSEA")) {
  dir.create(
    "Figures/GSEA",
    recursive = TRUE
  )
}

png(
  "Figures/GSEA/Hallmark_GSEA_Top_Pathways.png",
  width = 2600,
  height = 2000,
  res = 300
)

par(
  mar = c(5, 12, 4, 2)
)

plot(
  gsea_plot_data$NES,
  seq_len(nrow(gsea_plot_data)),
  type = "n",
  yaxt = "n",
  xlab = "Normalized Enrichment Score (NES)",
  ylab = "",
  main = "Hallmark GSEA: TNBC vs Normal",
  xlim = range(
    gsea_plot_data$NES
  ) + c(-0.3, 0.3)
)

axis(
  2,
  at = seq_len(nrow(gsea_plot_data)),
  labels = gsea_plot_data$Pathway,
  las = 1,
  cex.axis = 0.75
)

abline(
  v = 0,
  lty = 2
)

segments(
  x0 = 0,
  y0 = seq_len(nrow(gsea_plot_data)),
  x1 = gsea_plot_data$NES,
  y1 = seq_len(nrow(gsea_plot_data)),
  lwd = 5
)

points(
  gsea_plot_data$NES,
  seq_len(nrow(gsea_plot_data)),
  pch = 19,
  cex = 1.4
)

dev.off()

par(
  mfrow = c(1, 1)
)


# ============================================================
# 12. Final summary
# ============================================================

cat(
  "GSEA analysis completed.\n"
)

cat(
  "Ranked genes:",
  length(gsea_symbol_list),
  "\n"
)

cat(
  "Hallmark pathways:",
  length(hallmark_pathways),
  "\n"
)

cat(
  "Significant pathways:",
  nrow(gsea_significant),
  "\n"
)

cat(
  "TNBC-enriched pathways:",
  nrow(gsea_tnbc),
  "\n"
)

cat(
  "Normal-enriched pathways:",
  nrow(gsea_normal),
  "\n"
)

cat(
  "GSEA figure saved to:",
  "Figures/GSEA/Hallmark_GSEA_Top_Pathways.png",
  "\n"
)