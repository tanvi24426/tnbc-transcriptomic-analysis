# ============================================================
# 10_gene_pathway_integration.R
# Candidate Gene–Pathway Integration
# ============================================================


# ============================================================
# 1. Load required packages
# ============================================================

library(pheatmap)


# ============================================================
# 2. Load final candidate-gene panel
# ============================================================

final_genes <- read.csv(
  "Results/Final_10_Gene_Shortlist.csv",
  stringsAsFactors = FALSE
)

candidate_symbols <- final_genes$Symbol


# ============================================================
# 3. Load GSEA results
# ============================================================

gsea_results <- read.csv(
  "Results/GSEA/Hallmark_GSEA_All_Results.csv",
  stringsAsFactors = FALSE
)

# Keep significant Hallmark pathways.
significant_gsea <- gsea_results[
  gsea_results$padj < 0.05,
  ,
  drop = FALSE
]


# ============================================================
# 4. Recover leading-edge genes
# ============================================================

# The leadingEdge column was saved as a semicolon-separated
# character string in the CSV.

leading_edge_lists <- strsplit(
  significant_gsea$leadingEdge,
  split = ";",
  fixed = TRUE
)

names(leading_edge_lists) <- significant_gsea$pathway


# ============================================================
# 5. Calculate candidate-gene/pathway overlap
# ============================================================

candidate_pathway_overlap <- lapply(
  leading_edge_lists,
  function(leading_genes) {
    
    intersect(
      candidate_symbols,
      leading_genes
    )
    
  }
)

# Keep pathways containing at least one candidate gene.
candidate_pathway_overlap <- candidate_pathway_overlap[
  lengths(candidate_pathway_overlap) > 0
]


# ============================================================
# 6. Create long-format integration table
# ============================================================

integration_rows <- list()

row_counter <- 1

for (pathway in names(candidate_pathway_overlap)) {
  
  genes <- candidate_pathway_overlap[[pathway]]
  
  for (gene in genes) {
    
    pathway_info <- significant_gsea[
      significant_gsea$pathway == pathway,
      ,
      drop = FALSE
    ]
    
    gene_info <- final_genes[
      final_genes$Symbol == gene,
      ,
      drop = FALSE
    ]
    
    integration_rows[[row_counter]] <- data.frame(
      Gene = gene,
      Pathway = pathway,
      NES = pathway_info$NES[1],
      padj = pathway_info$padj[1],
      log2FoldChange = gene_info$log2FoldChange[1],
      Direction = ifelse(
        pathway_info$NES[1] > 0,
        "TNBC-enriched",
        "Normal-enriched"
      ),
      stringsAsFactors = FALSE
    )
    
    row_counter <- row_counter + 1
    
  }
}

gene_pathway_integration <- do.call(
  rbind,
  integration_rows
)


# ============================================================
# 7. Add pathway overlap counts
# ============================================================

pathway_overlap_counts <- aggregate(
  Gene ~ Pathway,
  data = gene_pathway_integration,
  FUN = length
)

colnames(pathway_overlap_counts)[2] <-
  "Candidate_Gene_Count"

gene_pathway_integration <- merge(
  gene_pathway_integration,
  pathway_overlap_counts,
  by = "Pathway",
  all.x = TRUE
)

# Restore a useful ordering.
gene_pathway_integration <- gene_pathway_integration[
  order(
    gene_pathway_integration$NES,
    gene_pathway_integration$padj,
    gene_pathway_integration$Gene
  ),
  ,
  drop = FALSE
]


# ============================================================
# 8. Save integration results
# ============================================================

if (!dir.exists("Results/Gene_Pathway_Integration")) {
  
  dir.create(
    "Results/Gene_Pathway_Integration",
    recursive = TRUE
  )
  
}

write.csv(
  gene_pathway_integration,
  "Results/Gene_Pathway_Integration/Candidate_Gene_Pathway_Overlap.csv",
  row.names = FALSE
)


# ============================================================
# 9. Create candidate × pathway matrix
# ============================================================

# Initialize a matrix of zeros.
integration_matrix <- matrix(
  0,
  nrow = length(candidate_symbols),
  ncol = length(candidate_pathway_overlap)
)

rownames(integration_matrix) <- candidate_symbols

colnames(integration_matrix) <-
  names(candidate_pathway_overlap)

# Mark candidate genes that occur in pathway leading edges.
for (pathway in names(candidate_pathway_overlap)) {
  
  genes <- candidate_pathway_overlap[[pathway]]
  
  integration_matrix[
    intersect(
      candidate_symbols,
      genes
    ),
    pathway
  ] <- 1
  
}


# ============================================================
# 10. Order genes by biological program
# ============================================================

program_order <- c(
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

program_order <- intersect(
  program_order,
  rownames(integration_matrix)
)

integration_matrix <- integration_matrix[
  program_order,
  ,
  drop = FALSE
]


# ============================================================
# 11. Order pathways by NES
# ============================================================

pathway_order <- significant_gsea$pathway[
  significant_gsea$pathway %in%
    colnames(integration_matrix)
]

pathway_order <- pathway_order[
  order(
    significant_gsea$NES[
      match(
        pathway_order,
        significant_gsea$pathway
      )
    ]
  )
]

integration_matrix <- integration_matrix[
  ,
  pathway_order,
  drop = FALSE
]


# ============================================================
# 12. Create heatmap
# ============================================================

if (!dir.exists("Figures/Gene_Pathway_Integration")) {
  
  dir.create(
    "Figures/Gene_Pathway_Integration",
    recursive = TRUE
  )
  
}

png(
  "Figures/Gene_Pathway_Integration/Candidate_Gene_Pathway_Heatmap.png",
  width = 3600,
  height = 2400,
  res = 300
)

# Create cleaner pathway labels.
clean_pathway_names <- gsub(
  "^HALLMARK_",
  "",
  colnames(integration_matrix)
)

clean_pathway_names <- gsub(
  "_",
  " ",
  clean_pathway_names
)

colnames(integration_matrix) <- clean_pathway_names

pheatmap(
  integration_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  fontsize = 10,
  fontsize_row = 11,
  fontsize_col = 8,
  angle_col = 45,
  main = "Candidate Genes in GSEA Leading-Edge Pathways",
  border_color = NA
)

dev.off()


# ============================================================
# 13. Summary statistics
# ============================================================

cat(
  "Candidate gene–pathway integration completed.\n"
)

cat(
  "Candidate genes:",
  length(candidate_symbols),
  "\n"
)

cat(
  "Significant GSEA pathways:",
  nrow(significant_gsea),
  "\n"
)

cat(
  "Pathways containing candidate genes:",
  length(candidate_pathway_overlap),
  "\n"
)

cat(
  "Candidate-pathway associations:",
  nrow(gene_pathway_integration),
  "\n"
)

cat(
  "Integration results saved to:",
  "Results/Gene_Pathway_Integration/Candidate_Gene_Pathway_Overlap.csv",
  "\n"
)

cat(
  "Heatmap saved to:",
  "Figures/Gene_Pathway_Integration/Candidate_Gene_Pathway_Heatmap.png",
  "\n"
)