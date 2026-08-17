# ============================================================
# 04_candidate_prioritization.R
# Candidate Gene Prioritization
# ============================================================


# ------------------------------------------------------------
# 1. Define stringent candidate-gene criteria
# ------------------------------------------------------------

# Candidate genes were selected using:
# adjusted p-value < 0.05
# absolute log2 fold change >= 2
# baseMean >= 50
#
# These criteria identify genes with strong statistical evidence,
# substantial expression, and a large expression change.

candidate_genes <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 2 &
    res_df$baseMean >= 50,
]

# Display the number of candidate genes
nrow(candidate_genes)


# ------------------------------------------------------------
# 2. Separate upregulated and downregulated candidates
# ------------------------------------------------------------

candidate_up <- candidate_genes[
  candidate_genes$log2FoldChange >= 2,
]

candidate_down <- candidate_genes[
  candidate_genes$log2FoldChange <= -2,
]

# Count candidates in each direction
nrow(candidate_up)
nrow(candidate_down)


# ------------------------------------------------------------
# 3. Annotate candidate genes
# ------------------------------------------------------------

# Add gene symbols and functional annotations
# using the human annotation table

candidate_genes_annotated <- merge(
  candidate_genes,
  annotation[
    ,
    c(
      "GeneID",
      "Symbol",
      "Description",
      "GeneType",
      "EnsemblGeneID"
    )
  ],
  by = "GeneID",
  all.x = TRUE
)

# Sort candidates by adjusted p-value
candidate_genes_annotated <- candidate_genes_annotated[
  order(candidate_genes_annotated$padj),
]

# Display the first few annotated candidates
head(candidate_genes_annotated)


# ------------------------------------------------------------
# 4. Rank upregulated candidates
# ------------------------------------------------------------

# Sort upregulated candidates by adjusted p-value
candidate_up_sorted <- candidate_up[
  order(candidate_up$padj),
]

# Add annotation
candidate_up_annotated <- merge(
  candidate_up_sorted,
  annotation[
    ,
    c(
      "GeneID",
      "Symbol",
      "Description",
      "GeneType",
      "EnsemblGeneID"
    )
  ],
  by = "GeneID",
  all.x = TRUE
)

# Sort again after merging
candidate_up_annotated <- candidate_up_annotated[
  order(candidate_up_annotated$padj),
]

# Display the strongest upregulated candidates
head(
  candidate_up_annotated,
  20
)


# ------------------------------------------------------------
# 5. Rank downregulated candidates
# ------------------------------------------------------------

# Sort downregulated candidates by adjusted p-value
candidate_down_sorted <- candidate_down[
  order(candidate_down$padj),
]

# Add annotation
candidate_down_annotated <- merge(
  candidate_down_sorted,
  annotation[
    ,
    c(
      "GeneID",
      "Symbol",
      "Description",
      "GeneType",
      "EnsemblGeneID"
    )
  ],
  by = "GeneID",
  all.x = TRUE
)

# Sort again after merging
candidate_down_annotated <- candidate_down_annotated[
  order(candidate_down_annotated$padj),
]

# Display the strongest downregulated candidates
head(
  candidate_down_annotated,
  20
)


# ------------------------------------------------------------
# 6. Identify pathway-associated candidates
# ------------------------------------------------------------

# Extract genes represented in the enriched KEGG pathways
# from the upregulated and downregulated enrichment results

kegg_up_genes <- unique(
  unlist(
    strsplit(
      kegg_up_df$geneID,
      "/"
    )
  )
)

kegg_down_genes <- unique(
  unlist(
    strsplit(
      kegg_down_df$geneID,
      "/"
    )
  )
)

# Identify upregulated candidates that occur in
# enriched KEGG pathways
candidate_up_pathway <- candidate_up_annotated[
  candidate_up_annotated$GeneID %in% kegg_up_genes,
]

# Identify downregulated candidates that occur in
# enriched KEGG pathways
candidate_down_pathway <- candidate_down_annotated[
  candidate_down_annotated$GeneID %in% kegg_down_genes,
]

# Sort both groups by adjusted p-value
candidate_up_pathway <- candidate_up_pathway[
  order(candidate_up_pathway$padj),
]

candidate_down_pathway <- candidate_down_pathway[
  order(candidate_down_pathway$padj),
]

# Display the strongest pathway-associated candidates
head(
  candidate_up_pathway,
  20
)

head(
  candidate_down_pathway,
  20
)


# ------------------------------------------------------------
# 7. Create review tables
# ------------------------------------------------------------

# Select the most informative columns for reviewing
# the strongest upregulated candidates

up_shortlist_review <- candidate_up_pathway[
  ,
  c(
    "Symbol",
    "Description",
    "baseMean",
    "log2FoldChange",
    "padj",
    "GeneType"
  )
]

# Select the most informative columns for reviewing
# the strongest downregulated candidates

down_shortlist_review <- candidate_down_pathway[
  ,
  c(
    "Symbol",
    "Description",
    "baseMean",
    "log2FoldChange",
    "padj",
    "GeneType"
  )
]

# Display the top candidates for review
head(
  up_shortlist_review,
  20
)

head(
  down_shortlist_review,
  20
)


# ------------------------------------------------------------
# 8. Define the final 10 prioritized candidates
# ------------------------------------------------------------

# Five upregulated genes were selected to represent the
# strongest cell-cycle/proliferative expression program.

final_up_symbols <- c(
  "CDC20",
  "BUB1",
  "TRIP13",
  "PLK1",
  "AURKB"
)

# Five downregulated genes were selected to represent the
# strongest lipid/metabolic expression program.

final_down_symbols <- c(
  "PNPLA2",
  "PPARG",
  "LIPE",
  "LEP",
  "CIDEC"
)

# Extract the final candidates from the annotated tables

final_up <- candidate_up_annotated[
  candidate_up_annotated$Symbol %in% final_up_symbols,
]

final_down <- candidate_down_annotated[
  candidate_down_annotated$Symbol %in% final_down_symbols,
]


# ------------------------------------------------------------
# 9. Combine the final candidates
# ------------------------------------------------------------

final_candidates <- rbind(
  final_up,
  final_down
)

# Order by direction and statistical significance
final_candidates <- final_candidates[
  order(
    -(
      final_candidates$log2FoldChange > 0
    ),
    final_candidates$padj
  ),
]

# Display the final candidate list
final_candidates[
  ,
  c(
    "Symbol",
    "Description",
    "baseMean",
    "log2FoldChange",
    "padj",
    "GeneType",
    "EnsemblGeneID"
  )
]


# ------------------------------------------------------------
# 10. Create a simplified final results table
# ------------------------------------------------------------

final_results_table <- final_candidates[
  ,
  c(
    "Symbol",
    "log2FoldChange",
    "padj",
    "baseMean",
    "Description"
  )
]

# Add expression direction
final_results_table$Direction <- ifelse(
  final_results_table$log2FoldChange > 0,
  "Upregulated",
  "Downregulated"
)

# Reorder columns
final_results_table <- final_results_table[
  ,
  c(
    "Symbol",
    "Direction",
    "log2FoldChange",
    "padj",
    "baseMean",
    "Description"
  )
]

# Display the final results table
final_results_table


# ------------------------------------------------------------
# 11. Save candidate-gene results
# ------------------------------------------------------------

# Create Results folder if necessary
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save all annotated candidate genes
write.csv(
  candidate_genes_annotated,
  "Results/All_Candidate_Genes.csv",
  row.names = FALSE
)

# Save pathway-associated candidates
write.csv(
  candidate_up_pathway,
  "Results/Pathway_Associated_Upregulated_Candidates.csv",
  row.names = FALSE
)

write.csv(
  candidate_down_pathway,
  "Results/Pathway_Associated_Downregulated_Candidates.csv",
  row.names = FALSE
)

# Save the final 10-gene shortlist
write.csv(
  final_candidates,
  "Results/Final_10_Gene_Shortlist.csv",
  row.names = FALSE
)

# Save the simplified final results table
write.csv(
  final_results_table,
  "Results/Final_10_Gene_Results_Table.csv",
  row.names = FALSE
)
