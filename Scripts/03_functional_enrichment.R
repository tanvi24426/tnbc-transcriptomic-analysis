# ============================================================
# 03_functional_enrichment.R
# Functional enrichment analysis: GO and KEGG
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)


# ------------------------------------------------------------
# 2. Prepare gene lists for enrichment analysis
# ------------------------------------------------------------

# Extract GeneIDs for upregulated genes
up_ids <- unique(
  as.character(upregulated$GeneID)
)

# Extract GeneIDs for downregulated genes
down_ids <- unique(
  as.character(downregulated$GeneID)
)

# Extract all genes tested in the differential expression analysis
# This will be used as the background/universe
background_ids <- unique(
  as.character(res_df$GeneID)
)


# ------------------------------------------------------------
# 3. GO Biological Process enrichment
# ------------------------------------------------------------

# Perform GO Biological Process enrichment
# for the upregulated genes
go_up <- enrichGO(
  gene = up_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Perform GO Biological Process enrichment
# for the downregulated genes
go_down <- enrichGO(
  gene = down_ids,
  universe = background_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)


# ------------------------------------------------------------
# 4. Inspect GO enrichment results
# ------------------------------------------------------------

# Convert GO results to regular data frames
go_up_df <- as.data.frame(go_up)
go_down_df <- as.data.frame(go_down)

# Display the top enriched biological processes
head(go_up_df, 10)
head(go_down_df, 10)


# ------------------------------------------------------------
# 5. KEGG pathway enrichment
# ------------------------------------------------------------

# Perform KEGG pathway enrichment for upregulated genes
kegg_up <- enrichKEGG(
  gene = up_ids,
  universe = background_ids,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Perform KEGG pathway enrichment for downregulated genes
kegg_down <- enrichKEGG(
  gene = down_ids,
  universe = background_ids,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)


# ------------------------------------------------------------
# 6. Inspect KEGG enrichment results
# ------------------------------------------------------------

# Convert KEGG results to regular data frames
kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

# Display the top 10 enriched KEGG pathways
head(kegg_up_df, 10)
head(kegg_down_df, 10)


# ------------------------------------------------------------
# 7. Save enrichment results
# ------------------------------------------------------------

# Create Results folder if it does not already exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save GO enrichment results
write.csv(
  go_up_df,
  "Results/GO_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  go_down_df,
  "Results/GO_Downregulated.csv",
  row.names = FALSE
)

# Save KEGG enrichment results
write.csv(
  kegg_up_df,
  "Results/KEGG_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  kegg_down_df,
  "Results/KEGG_Downregulated.csv",
  row.names = FALSE
)