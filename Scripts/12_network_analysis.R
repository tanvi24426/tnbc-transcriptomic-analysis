# ============================================================
# 12_network_analysis.R
# STRING PPI Network Analysis of Final Candidate Genes
# ============================================================


# ============================================================
# 1. Load required package
# ============================================================

library(igraph)


# ============================================================
# 2. Load final candidate genes
# ============================================================

final_genes <- read.csv(
  "Results/Final_10_Gene_Shortlist.csv",
  stringsAsFactors = FALSE
)

candidate_symbols <- unique(
  final_genes$Symbol
)

print(candidate_symbols)


# ============================================================
# 3. Define STRING API settings
# ============================================================

string_base_url <- "https://string-db.org/api"

species <- 9606

required_score <- 700

caller_identity <- "tnbc-transcriptomic-analysis"


# ============================================================
# 4. Map gene symbols to STRING identifiers
# ============================================================

identifier_string <- paste(
  candidate_symbols,
  collapse = "%0d"
)

mapping_url <- paste0(
  string_base_url,
  "/tsv/get_string_ids?",
  "identifiers=",
  identifier_string,
  "&species=",
  species,
  "&echo_query=1",
  "&caller_identity=",
  URLencode(
    caller_identity,
    reserved = TRUE
  )
)

mapping_file <- tempfile(
  fileext = ".tsv"
)

download.file(
  mapping_url,
  mapping_file,
  mode = "wb",
  quiet = TRUE
)

string_mapping <- read.delim(
  mapping_file,
  stringsAsFactors = FALSE
)

print(string_mapping)


# ============================================================
# 5. Check mapping completeness
# ============================================================

mapped_genes <- unique(
  string_mapping$preferredName
)

unmapped_genes <- setdiff(
  candidate_symbols,
  mapped_genes
)

cat(
  "Candidate genes:",
  length(candidate_symbols),
  "\n"
)

cat(
  "Mapped genes:",
  length(mapped_genes),
  "\n"
)

cat(
  "Unmapped genes:",
  length(unmapped_genes),
  "\n"
)

if (
  length(unmapped_genes) > 0
) {
  
  cat(
    "Unmapped:",
    paste(
      unmapped_genes,
      collapse = ", "
    ),
    "\n"
  )
  
}


# ============================================================
# 6. Prepare STRING identifiers
# ============================================================

string_ids <- unique(
  string_mapping$stringId
)

network_identifier_string <- paste(
  string_ids,
  collapse = "%0d"
)


# ============================================================
# 7. Retrieve STRING interaction network
# ============================================================

network_url <- paste0(
  string_base_url,
  "/tsv/network?",
  "identifiers=",
  network_identifier_string,
  "&species=",
  species,
  "&required_score=",
  required_score,
  "&network_type=functional",
  "&caller_identity=",
  URLencode(
    caller_identity,
    reserved = TRUE
  )
)

network_file <- tempfile(
  fileext = ".tsv"
)

download.file(
  network_url,
  network_file,
  mode = "wb",
  quiet = TRUE
)

ppi_network <- read.delim(
  network_file,
  stringsAsFactors = FALSE
)

cat(
  "STRING interactions retrieved:",
  nrow(ppi_network),
  "\n"
)

print(
  head(ppi_network)
)


# ============================================================
# 8. Restrict interactions to the 10 candidate proteins
# ============================================================

candidate_string_ids <- unique(
  string_mapping$stringId
)

candidate_network <- ppi_network[
  ppi_network$stringId_A %in% candidate_string_ids &
    ppi_network$stringId_B %in% candidate_string_ids,
  ,
  drop = FALSE
]


# ============================================================
# 9. Replace STRING IDs with gene symbols
# ============================================================

id_to_gene <- setNames(
  string_mapping$preferredName,
  string_mapping$stringId
)

candidate_network$preferredName_A <- unname(
  id_to_gene[
    candidate_network$stringId_A
  ]
)

candidate_network$preferredName_B <- unname(
  id_to_gene[
    candidate_network$stringId_B
  ]
)

cat(
  "High-confidence candidate interactions:",
  nrow(candidate_network),
  "\n"
)

print(
  candidate_network[
    ,
    c(
      "preferredName_A",
      "preferredName_B",
      "score"
    )
  ]
)


# ============================================================
# 10. Save raw candidate network
# ============================================================

if (
  !dir.exists("Results/Network")
) {
  
  dir.create(
    "Results/Network",
    recursive = TRUE
  )
  
}

if (
  !dir.exists("Figures/Network")
) {
  
  dir.create(
    "Figures/Network",
    recursive = TRUE
  )
  
}

write.csv(
  candidate_network,
  "Results/Network/Candidate_PPI_Network.csv",
  row.names = FALSE
)


# ============================================================
# 11. Create igraph network
# ============================================================

if (
  nrow(candidate_network) > 0
) {
  
  edge_table <- candidate_network[
    ,
    c(
      "preferredName_A",
      "preferredName_B",
      "score"
    ),
    drop = FALSE
  ]
  
  colnames(edge_table) <- c(
    "from",
    "to",
    "score"
  )
  
  graph <- graph_from_data_frame(
    edge_table,
    directed = FALSE
  )
  
} else {
  
  graph <- make_empty_graph()
  
}


# ============================================================
# 12. Calculate network centrality
# ============================================================

if (
  vcount(graph) > 0
) {
  
  degree_values <- degree(
    graph,
    mode = "all"
  )
  
  betweenness_values <- betweenness(
    graph,
    directed = FALSE,
    normalized = TRUE
  )
  
  closeness_values <- closeness(
    graph,
    normalized = TRUE
  )
  
  centrality_results <- data.frame(
    
    Gene = names(
      degree_values
    ),
    
    Degree = as.numeric(
      degree_values
    ),
    
    Betweenness = as.numeric(
      betweenness_values
    ),
    
    Closeness = as.numeric(
      closeness_values
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  centrality_results <- centrality_results[
    order(
      -centrality_results$Degree,
      -centrality_results$Betweenness
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  centrality_results <- data.frame(
    Gene = candidate_symbols,
    Degree = 0,
    Betweenness = 0,
    Closeness = 0,
    stringsAsFactors = FALSE
  )
  
}


# ============================================================
# 13. Add biological program annotation
# ============================================================

centrality_results$Program <- ifelse(
  
  centrality_results$Gene %in%
    c(
      "CDC20",
      "BUB1",
      "TRIP13",
      "PLK1",
      "AURKB"
    ),
  
  "Proliferation",
  
  ifelse(
    
    centrality_results$Gene %in%
      c(
        "PNPLA2",
        "PPARG",
        "LIPE",
        "LEP",
        "CIDEC"
      ),
    
    "Metabolic",
    
    "Other"
    
  )
  
)


# ============================================================
# 14. Save centrality results
# ============================================================

write.csv(
  centrality_results,
  "Results/Network/Candidate_Network_Centrality.csv",
  row.names = FALSE
)


# ============================================================
# 15. STRING PPI enrichment
# ============================================================

ppi_enrichment_url <- paste0(
  string_base_url,
  "/tsv/ppi_enrichment?",
  "identifiers=",
  network_identifier_string,
  "&species=",
  species,
  "&required_score=",
  required_score,
  "&caller_identity=",
  URLencode(
    caller_identity,
    reserved = TRUE
  )
)

enrichment_file <- tempfile(
  fileext = ".tsv"
)

download.file(
  ppi_enrichment_url,
  enrichment_file,
  mode = "wb",
  quiet = TRUE
)

ppi_enrichment <- read.delim(
  enrichment_file,
  stringsAsFactors = FALSE
)

write.csv(
  ppi_enrichment,
  "Results/Network/Candidate_PPI_Enrichment.csv",
  row.names = FALSE
)

# ============================================================
# 16. Improved PPI Network Figure
# ============================================================

png(
  "Figures/Network/Candidate_PPI_Network.png",
  width = 3600,
  height = 3000,
  res = 300
)

# ------------------------------------------------------------
# Biological modules
# ------------------------------------------------------------

proliferation_genes <- c(
  "CDC20", "BUB1", "TRIP13", "PLK1", "AURKB"
)

metabolic_genes <- c(
  "PNPLA2", "PPARG", "LIPE", "LEP", "CIDEC"
)

# ------------------------------------------------------------
# Manual layout
# ------------------------------------------------------------

layout_matrix <- matrix(
  c(
    # Cell-cycle / proliferation
    -1.00,  0.00,     # TRIP13
    -0.45,  0.60,     # CDC20
    0.20,  0.48,     # PLK1
    0.25, -0.38,     # AURKB
    -0.45, -0.58,     # BUB1
    
    # Metabolic / adipose
    0.85,  0.62,     # PNPLA2
    1.30,  0.62,     # LIPE
    0.78,  0.05,     # PPARG
    1.50,  0.05,     # LEP
    1.20, -0.48      # CIDEC
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(layout_matrix) <- c(
  "TRIP13",
  "CDC20",
  "PLK1",
  "AURKB",
  "BUB1",
  "PNPLA2",
  "LIPE",
  "PPARG",
  "LEP",
  "CIDEC"
)

layout_matrix <- layout_matrix[
  match(
    V(graph)$name,
    rownames(layout_matrix)
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Network statistics
# ------------------------------------------------------------

degree_values <- degree(
  graph,
  mode = "all"
)

edge_scores <- E(graph)$score

# ------------------------------------------------------------
# Node colours
# ------------------------------------------------------------

node_colors <- ifelse(
  V(graph)$name %in% proliferation_genes,
  "indianred1",
  "skyblue2"
)

# ------------------------------------------------------------
# Node sizes
# ------------------------------------------------------------

node_sizes <- 18 + 4 * degree_values

# ------------------------------------------------------------
# Edge widths
# ------------------------------------------------------------

edge_widths <- 1.5 +
  5 * (
    (edge_scores - min(edge_scores)) /
      (max(edge_scores) - min(edge_scores))
  )

# ------------------------------------------------------------
# Create plotting area
# ------------------------------------------------------------

plot.new()

plot.window(
  xlim = c(-1.45, 1.95),
  ylim = c(-1.25, 1.10),
  asp = 1
)

# ------------------------------------------------------------
# Module backgrounds
# ------------------------------------------------------------

polygon(
  c(
    -1.35, -0.90, -0.25,
    0.40,  0.45,  0.15,
    -0.30, -0.90, -1.35
  ),
  c(
    0.00,  0.78,  0.83,
    0.40, -0.25, -0.72,
    -0.82, -0.70,  0.00
  ),
  col = "mistyrose",
  border = NA
)

polygon(
  c(
    0.55,  0.90,  1.45,
    1.78,  1.78,  1.35,
    0.72,  0.50
  ),
  c(
    0.05,  0.78,  0.78,
    0.35, -0.35, -0.72,
    -0.72, -0.20
  ),
  col = "lightcyan",
  border = NA
)

# ------------------------------------------------------------
# Module titles
# ------------------------------------------------------------

text(
  -0.45,
  0.96,
  "Cell-cycle / Proliferation Module",
  font = 2,
  cex = 1.15
)

text(
  1.15,
  0.96,
  "Metabolic / Adipose Module",
  font = 2,
  cex = 1.15
)

# ------------------------------------------------------------
# Draw the ACTUAL network
# ------------------------------------------------------------

plot(
  graph,
  
  layout = layout_matrix,
  
  rescale = FALSE,
  
  xlim = c(-1.45, 1.95),
  ylim = c(-1.25, 1.10),
  
  add = TRUE,
  
  vertex.color = node_colors,
  
  vertex.frame.color = "grey25",
  
  vertex.frame.width = 1,
  
  vertex.size = node_sizes,
  
  vertex.label = V(graph)$name,
  
  vertex.label.color = "black",
  
  vertex.label.cex = 1.05,
  
  vertex.label.family = "sans",
  
  edge.color = "grey45",
  
  edge.width = edge_widths,
  
  edge.label = NA,
  
  edge.curved = 0.04
)

# ------------------------------------------------------------
# Add edge confidence scores manually
#
# This prevents the numbers from being placed directly
# on top of the network lines.
# ------------------------------------------------------------

edge_ends <- ends(
  graph,
  E(graph),
  names = TRUE
)

for (i in seq_len(nrow(edge_ends))) {
  
  gene1 <- edge_ends[i, 1]
  gene2 <- edge_ends[i, 2]
  
  p1 <- layout_matrix[gene1, ]
  p2 <- layout_matrix[gene2, ]
  
  # midpoint
  mid_x <- (p1[1] + p2[1]) / 2
  mid_y <- (p1[2] + p2[2]) / 2
  
  # perpendicular offset
  dx <- p2[1] - p1[1]
  dy <- p2[2] - p1[2]
  
  length_edge <- sqrt(dx^2 + dy^2)
  
  if (length_edge > 0) {
    
    offset_x <- -dy / length_edge * 0.045
    offset_y <-  dx / length_edge * 0.045
    
    text(
      mid_x + offset_x,
      mid_y + offset_y,
      labels = sprintf(
        "%.3f",
        edge_scores[i]
      ),
      cex = 0.72,
      col = "black",
      font = 1
    )
  }
}

# ------------------------------------------------------------
# Main title
# ------------------------------------------------------------

title(
  main =
    "Protein–Protein Interaction (PPI) Network of 10 Candidate Genes",
  font.main = 2,
  cex.main = 1.35
)

# ------------------------------------------------------------
# Subtitle
# ------------------------------------------------------------

mtext(
  "STRING high-confidence functional associations (combined score ≥ 700)",
  side = 3,
  line = 0.2,
  cex = 1.0,
  font = 2
)

# ------------------------------------------------------------
# STRING confidence legend
# ------------------------------------------------------------

legend(
  "topright",
  
  legend = c(
    "0.900–1.000",
    "0.800–0.899",
    "0.700–0.799"
  ),
  
  lwd = c(
    5,
    3.5,
    2
  ),
  
  col = "grey45",
  
  bty = "n",
  
  title = "STRING confidence",
  
  cex = 0.80
)

# ------------------------------------------------------------
# Biological module legend
# ------------------------------------------------------------

legend(
  "bottomleft",
  
  legend = c(
    "Cell-cycle / Proliferation",
    "Metabolic / Adipose"
  ),
  
  pch = 21,
  
  pt.bg = c(
    "indianred1",
    "skyblue2"
  ),
  
  pt.cex = 1.5,
  
  bty = "n",
  
  title = "Biological module",
  
  cex = 0.82
)

# ------------------------------------------------------------
# Explanation / interpretation box
# ------------------------------------------------------------

rect(
  -1.35,
  -1.20,
  1.78,
  -0.88,
  col = "white",
  border = "grey60",
  lwd = 1
)

# What the figure represents
text(
  0.20,
  -0.96,
  labels =
    "Nodes = prioritized candidate genes   |   Edges = STRING functional associations   |   Numbers = STRING confidence score",
  cex = 0.76
)

# Network summary
text(
  0.20,
  -1.04,
  labels = paste0(
    "Network summary: ",
    vcount(graph),
    " genes | ",
    ecount(graph),
    " high-confidence interactions | ",
    components(graph)$no,
    " connected components"
  ),
  cex = 0.76
)

# Biological interpretation
text(
  0.20,
  -1.12,
  labels =
    "Two major functional modules are observed: cell-cycle/proliferation and metabolic/adipose-related processes.",
  cex = 0.76
)

# ------------------------------------------------------------
# Close figure
# ------------------------------------------------------------

dev.off()

# ============================================================
# 17. Summary
# ============================================================

cat(
  "\nNetwork analysis completed.\n"
)

cat(
  "Candidate genes:",
  length(candidate_symbols),
  "\n"
)

cat(
  "Mapped STRING proteins:",
  length(string_ids),
  "\n"
)

cat(
  "High-confidence candidate interactions:",
  nrow(candidate_network),
  "\n"
)

cat(
  "Network nodes:",
  vcount(graph),
  "\n"
)

cat(
  "Network edges:",
  ecount(graph),
  "\n"
)

cat(
  "Results saved to:",
  "Results/Network/",
  "\n"
)

cat(
  "Figure saved to:",
  "Figures/Network/Candidate_PPI_Network.png",
  "\n"
)