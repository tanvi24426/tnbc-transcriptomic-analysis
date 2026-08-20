# ============================================================
# 14. Drug-Target Analysis
# TNBC Transcriptomic Analysis
# ============================================================


# ============================================================
# Libraries
# ============================================================

library(httr2)
library(jsonlite)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# ============================================================
# Directories
# ============================================================

dir.create(
  "Results/Drug_Target",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/Drug_Target",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Candidate genes
# ============================================================

candidate_genes <- c(
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


# ============================================================
# Molecular program assignment
# ============================================================

gene_program <- data.frame(
  Gene = candidate_genes,
  Program = c(
    rep("Proliferation", 5),
    rep("Metabolic", 5)
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# DGIdb GraphQL API
# ============================================================

dgidb_url <- "https://dgidb.org/api/graphql"


# ============================================================
# Build GraphQL query
# ============================================================

gene_string <- paste(
  paste0('"', candidate_genes, '"'),
  collapse = ", "
)

graphql_query <- paste0(
  '{
    genes(names: [',
  gene_string,
  ']) {
      nodes {
        name
        conceptId
        interactions {
          drug {
            name
            conceptId
          }
          interactionScore
          interactionTypes {
            type
            directionality
          }
          interactionAttributes {
            name
            value
          }
          publications {
            pmid
          }
          sources {
            sourceDbName
          }
        }
      }
    }
  }'
)


# ============================================================
# Query DGIdb
# ============================================================

cat("\nQuerying DGIdb...\n")

request <- request(dgidb_url) |>
  req_method("POST") |>
  req_headers(
    `Content-Type` = "application/json"
  ) |>
  req_body_json(
    list(
      query = graphql_query
    )
  )

response <- req_perform(request)

response_status <- resp_status(response)

cat(
  "DGIdb HTTP status:",
  response_status,
  "\n"
)

if (response_status != 200) {
  stop(
    "DGIdb API request failed with HTTP status ",
    response_status
  )
}

response_data <- resp_body_json(
  response,
  simplifyVector = FALSE
)


# ============================================================
# Check API errors
# ============================================================

if (!is.null(response_data$errors)) {
  
  print(response_data$errors)
  
  stop(
    "DGIdb returned an API error."
  )
}


# ============================================================
# Extract gene nodes
# ============================================================

gene_nodes <- response_data$data$genes$nodes

cat(
  "Genes returned by DGIdb:",
  length(gene_nodes),
  "\n"
)


# ============================================================
# Extract drug-gene interactions
# ============================================================

interaction_list <- list()

counter <- 1

for (gene_node in gene_nodes) {
  
  if (is.null(gene_node$interactions)) {
    next
  }
  
  for (interaction in gene_node$interactions) {
    
    drug_name <- NA_character_
    drug_id <- NA_character_
    
    if (!is.null(interaction$drug)) {
      
      drug_name <- interaction$drug$name
      drug_id <- interaction$drug$conceptId
    }
    
    interaction_type <- NA_character_
    directionality <- NA_character_
    
    if (
      !is.null(interaction$interactionTypes) &&
      length(interaction$interactionTypes) > 0
    ) {
      
      interaction_type <- paste(
        vapply(
          interaction$interactionTypes,
          function(x) {
            
            if (is.null(x$type)) {
              NA_character_
            } else {
              x$type
            }
            
          },
          character(1)
        ),
        collapse = "; "
      )
      
      directionality <- paste(
        vapply(
          interaction$interactionTypes,
          function(x) {
            
            if (is.null(x$directionality)) {
              NA_character_
            } else {
              x$directionality
            }
            
          },
          character(1)
        ),
        collapse = "; "
      )
    }
    
    
    source_names <- NA_character_
    
    if (
      !is.null(interaction$sources) &&
      length(interaction$sources) > 0
    ) {
      
      source_names <- paste(
        unique(
          vapply(
            interaction$sources,
            function(x) {
              
              if (is.null(x$sourceDbName)) {
                NA_character_
              } else {
                x$sourceDbName
              }
              
            },
            character(1)
          )
        ),
        collapse = "; "
      )
    }
    
    
    pmids <- NA_character_
    
    if (
      !is.null(interaction$publications) &&
      length(interaction$publications) > 0
    ) {
      
      pmids <- paste(
        unique(
          vapply(
            interaction$publications,
            function(x) {
              
              if (is.null(x$pmid)) {
                NA_character_
              } else {
                as.character(x$pmid)
              }
              
            },
            character(1)
          )
        ),
        collapse = "; "
      )
    }
    
    
    interaction_list[[counter]] <- data.frame(
      Gene = gene_node$name,
      Drug = drug_name,
      Drug_Concept_ID = drug_id,
      Interaction_Score = interaction$interactionScore,
      Interaction_Type = interaction_type,
      Directionality = directionality,
      Sources = source_names,
      PubMed_IDs = pmids,
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
}


# ============================================================
# Combine interactions
# ============================================================

if (length(interaction_list) == 0) {
  
  stop(
    "No drug-gene interactions were returned by DGIdb."
  )
}

drug_target_results <- bind_rows(
  interaction_list
)


# ============================================================
# Remove incomplete records
# ============================================================

drug_target_results <- drug_target_results |>
  filter(
    !is.na(Gene),
    !is.na(Drug),
    Drug != ""
  )


# ============================================================
# Add molecular program
# ============================================================

drug_target_results <- drug_target_results |>
  left_join(
    gene_program,
    by = "Gene"
  )


# ============================================================
# Remove duplicate interactions
# ============================================================

drug_target_results <- drug_target_results |>
  distinct(
    Gene,
    Drug,
    Interaction_Type,
    Directionality,
    .keep_all = TRUE
  )


# ============================================================
# Summary
# ============================================================

cat(
  "\nDrug-target analysis completed.\n"
)

cat(
  "Unique candidate genes with interactions:",
  length(
    unique(drug_target_results$Gene)
  ),
  "\n"
)

cat(
  "Unique drugs identified:",
  length(
    unique(drug_target_results$Drug)
  ),
  "\n"
)

cat(
  "Drug-gene interaction records:",
  nrow(drug_target_results),
  "\n"
)


# ============================================================
# Gene-level summary
# ============================================================

gene_summary <- drug_target_results |>
  group_by(
    Gene,
    Program
  ) |>
  summarise(
    Number_of_Drugs = n_distinct(Drug),
    
    Number_of_Interactions = n(),
    
    Mean_Interaction_Score =
      mean(
        Interaction_Score,
        na.rm = TRUE
      ),
    
    Sources = paste(
      unique(
        unlist(
          strsplit(
            Sources,
            "; "
          )
        )
      ),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) |>
  arrange(
    desc(Number_of_Drugs)
  )


# ============================================================
# Drug-level summary
# ============================================================

drug_summary <- drug_target_results |>
  group_by(
    Drug
  ) |>
  summarise(
    
    Number_of_Target_Genes =
      n_distinct(Gene),
    
    Target_Genes =
      paste(
        sort(
          unique(Gene)
        ),
        collapse = ", "
      ),
    
    Programs =
      paste(
        sort(
          unique(Program)
        ),
        collapse = ", "
      ),
    
    Number_of_Sources =
      length(
        unique(
          unlist(
            strsplit(
              Sources,
              "; "
            )
          )
        )
      ),
    
    .groups = "drop"
    
  ) |>
  arrange(
    desc(Number_of_Target_Genes),
    Drug
  )


# ============================================================
# Save complete interaction table
# ============================================================

write.csv(
  drug_target_results,
  "Results/Drug_Target/Candidate_Drug_Targets.csv",
  row.names = FALSE
)


# ============================================================
# Save gene summary
# ============================================================

write.csv(
  gene_summary,
  "Results/Drug_Target/Candidate_Drug_Target_Summary.csv",
  row.names = FALSE
)


# ============================================================
# Save drug summary
# ============================================================

write.csv(
  drug_summary,
  "Results/Drug_Target/Drug_Target_Summary.csv",
  row.names = FALSE
)


# ============================================================
# Print gene summary
# ============================================================

cat(
  "\nTop genes by number of associated drugs:\n"
)

print(
  head(
    gene_summary,
    20
  )
)


# ============================================================
# Identify multi-target drugs
# ============================================================

multi_target_drugs <- drug_summary |>
  filter(
    Number_of_Target_Genes > 1
  )

cat(
  "\nDrugs associated with multiple candidate genes:\n"
)

print(
  multi_target_drugs
)


# ============================================================
# Save multi-target interaction table
# ============================================================

multi_target_summary <- drug_target_results |>
  filter(
    Drug %in% multi_target_drugs$Drug
  ) |>
  select(
    Gene,
    Program,
    Drug,
    Drug_Concept_ID,
    Interaction_Score,
    Interaction_Type,
    Directionality,
    Sources,
    PubMed_IDs
  ) |>
  arrange(
    Drug,
    Gene
  )

write.csv(
  multi_target_summary,
  "Results/Drug_Target/Multi_Target_Drug_Interactions.csv",
  row.names = FALSE
)


# ============================================================
# Gene-level drug association plot
# ============================================================

gene_summary_plot <- gene_summary |>
  arrange(
    Number_of_Drugs
  ) |>
  mutate(
    Gene = factor(
      Gene,
      levels = Gene
    )
  )

drug_count_plot <- ggplot(
  gene_summary_plot,
  aes(
    x = Gene,
    y = Number_of_Drugs,
    fill = Program
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  coord_flip() +
  
  labs(
    title = "DGIdb Drug Associations Across Candidate Genes",
    subtitle = "Number of recorded drug–gene associations by molecular program",
    x = "Candidate Gene",
    y = "Number of Associated Drugs",
    fill = "Molecular Program"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    legend.position = "bottom"
  )

print(
  drug_count_plot
)

ggsave(
  "Figures/Drug_Target/Candidate_Gene_Drug_Counts.png",
  drug_count_plot,
  width = 9,
  height = 6.5,
  dpi = 300
)


# ============================================================
# Focused multi-target drug network
# ============================================================

# Keep only Gene-Drug relationships
# Program information is already stored separately in
# gene_positions, so it is not included here.

multi_target_edges <- multi_target_summary |>
  select(
    Gene,
    Drug
  ) |>
  distinct()


# ============================================================
# Check multi-target network
# ============================================================

cat(
  "\nMulti-target drug-gene interactions:",
  nrow(multi_target_edges),
  "\n"
)

cat(
  "Multi-target drugs:",
  n_distinct(multi_target_edges$Drug),
  "\n"
)

cat(
  "Candidate genes represented:",
  n_distinct(multi_target_edges$Gene),
  "\n"
)


# ============================================================
# Define gene positions
# ============================================================

gene_positions <- data.frame(
  Gene = c(
    "PLK1",
    "AURKB",
    "PPARG",
    "LEP"
  ),
  
  x = 0,
  
  y = c(
    7,
    5.5,
    3.5,
    2
  ),
  
  Program = c(
    "Proliferation",
    "Proliferation",
    "Metabolic",
    "Metabolic"
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# Define drug positions
# ============================================================

drug_positions <- data.frame(
  Drug = c(
    "AMORFRUTIN A",
    "BENZBROMARONE",
    "CHEMBL:CHEMBL184450",
    "GW7647",
    "NVP-TAE684",
    "ESTRADIOL VALERATE",
    "TROGLITAZONE"
  ),
  
  x = 1,
  
  y = c(
    7.2,
    6.2,
    5.2,
    4.2,
    3.2,
    2.2,
    1.2
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# Build edge coordinates
# ============================================================

network_plot_data <- multi_target_edges |>
  
  left_join(
    gene_positions |>
      select(
        Gene,
        Gene_x = x,
        Gene_y = y
      ),
    by = "Gene"
  ) |>
  
  left_join(
    drug_positions |>
      select(
        Drug,
        Drug_x = x,
        Drug_y = y
      ),
    by = "Drug"
  )


# ============================================================
# Check for missing coordinates
# ============================================================

if (
  any(
    is.na(network_plot_data$Gene_x) |
    is.na(network_plot_data$Gene_y) |
    is.na(network_plot_data$Drug_x) |
    is.na(network_plot_data$Drug_y)
  )
) {
  
  warning(
    "Some drug-gene interactions do not have plotting coordinates."
  )
  
  network_plot_data <- network_plot_data |>
    filter(
      !is.na(Gene_x),
      !is.na(Gene_y),
      !is.na(Drug_x),
      !is.na(Drug_y)
    )
}


# ============================================================
# Prepare labels
# ============================================================

gene_plot_labels <- gene_positions |>
  mutate(
    Label = Gene
  )


drug_plot_labels <- drug_positions |>
  mutate(
    Label = str_wrap(
      Drug,
      width = 18
    )
  )


# ============================================================
# Focused network plot
# ============================================================

multi_target_network_plot <- ggplot() +
  
  # ----------------------------------------------------------
# Drug-gene interaction edges
# ----------------------------------------------------------

geom_segment(
  data = network_plot_data,
  aes(
    x = Gene_x,
    y = Gene_y,
    xend = Drug_x,
    yend = Drug_y
  ),
  linewidth = 1,
  alpha = 0.55,
  color = "grey50"
) +
  
  # ----------------------------------------------------------
# Candidate gene nodes
# ----------------------------------------------------------

geom_point(
  data = gene_positions,
  aes(
    x = x,
    y = y,
    fill = Program
  ),
  shape = 21,
  size = 10,
  stroke = 1.2,
  color = "black"
) +
  
  # ----------------------------------------------------------
# Drug nodes
# ----------------------------------------------------------

geom_point(
  data = drug_positions,
  aes(
    x = x,
    y = y
  ),
  shape = 21,
  size = 8,
  fill = "gold",
  stroke = 1.2,
  color = "black"
) +
  
  # ----------------------------------------------------------
# Gene labels
# ----------------------------------------------------------

geom_text(
  data = gene_plot_labels,
  aes(
    x = x - 0.06,
    y = y,
    label = Label
  ),
  hjust = 1,
  fontface = "bold",
  size = 5
) +
  
  # ----------------------------------------------------------
# Drug labels
# ----------------------------------------------------------

geom_text(
  data = drug_plot_labels,
  aes(
    x = x + 0.06,
    y = y,
    label = Label
  ),
  hjust = 0,
  size = 3.7,
  lineheight = 0.9
) +
  
  # ----------------------------------------------------------
# Column headings
# ----------------------------------------------------------

annotate(
  "text",
  x = 0,
  y = 8.0,
  label = "Candidate genes",
  fontface = "bold",
  size = 5
) +
  
  annotate(
    "text",
    x = 1,
    y = 8.0,
    label = "Multi-target drugs",
    fontface = "bold",
    size = 5
  ) +
  
  # ----------------------------------------------------------
# Program labels
# ----------------------------------------------------------

annotate(
  "text",
  x = -0.06,
  y = 6.25,
  label = "Proliferation",
  hjust = 1,
  size = 3.5,
  fontface = "italic"
) +
  
  annotate(
    "text",
    x = -0.06,
    y = 2.75,
    label = "Metabolic",
    hjust = 1,
    size = 3.5,
    fontface = "italic"
  ) +
  
  # ----------------------------------------------------------
# Titles
# ----------------------------------------------------------

labs(
  title = "Multi-Target Drug–Candidate Gene Network",
  
  subtitle =
    "DGIdb associations involving drugs linked to more than one prioritized candidate gene"
) +
  
  # ----------------------------------------------------------
# Gene program colors
# ----------------------------------------------------------

scale_fill_manual(
  values = c(
    "Proliferation" = "#D95F5F",
    "Metabolic" = "#4C8AC2"
  )
) +
  
  # ----------------------------------------------------------
# Plot limits
# ----------------------------------------------------------

scale_x_continuous(
  limits = c(
    -0.45,
    1.85
  ),
  expand = c(
    0,
    0
  )
) +
  
  scale_y_continuous(
    limits = c(
      0.5,
      8.4
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_void(
  base_size = 13
) +
  
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      color = "black"
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      margin = margin(b = 15),
      color = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    
    plot.margin = margin(
      20,
      30,
      20,
      30
    )
  )


# ============================================================
# Display network
# ============================================================

print(
  multi_target_network_plot
)


# ============================================================
# Save network
# ============================================================

ggsave(
  "Figures/Drug_Target/Multi_Target_Drug_Network.png",
  multi_target_network_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


# ============================================================
# Final message
# ============================================================

cat(
  "\nDrug/target analysis completed successfully.\n"
)

cat(
  "Results saved to: Results/Drug_Target/\n"
)

cat(
  "Figures saved to: Figures/Drug_Target/\n"
)

cat(
  "Multi-target network saved to: ",
  "Figures/Drug_Target/Multi_Target_Drug_Network.png\n"
)