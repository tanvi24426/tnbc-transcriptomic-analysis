# ============================================================
# 11_ROC_AUC_Analysis.R
# ROC / AUC Analysis of Candidate Genes and Molecular Signatures
# ============================================================


# ============================================================
# 1. Load required package
# ============================================================

library(pROC)


# ============================================================
# 2. Load final candidate genes
# ============================================================

final_genes <- read.csv(
  "Results/Final_10_Gene_Shortlist.csv",
  stringsAsFactors = FALSE
)

candidate_symbols <- final_genes$Symbol


# ============================================================
# 3. Load gene annotation
# ============================================================

annotation <- read.delim(
  gzfile("Data/Human.GRCh38.p13.annot.tsv.gz"),
  stringsAsFactors = FALSE
)

# Make sure GeneID is numeric.
annotation$GeneID <- as.integer(annotation$GeneID)


# ============================================================
# 4. Map candidate symbols to GeneIDs
# ============================================================

candidate_annotation <- annotation[
  annotation$Symbol %in% candidate_symbols,
  c("GeneID", "Symbol"),
  drop = FALSE
]

# Keep one record per candidate symbol.
candidate_annotation <- candidate_annotation[
  !duplicated(candidate_annotation$Symbol),
  ,
  drop = FALSE
]

candidate_annotation


# ============================================================
# 5. Load raw discovery counts
# ============================================================

counts <- read.delim(
  gzfile("Data/GSE233242_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 6. Load sample metadata
# ============================================================

sample_map <- read.csv(
  gzfile("Data/GSE233242_Sample_ID_to_GEO_ids.csv.gz"),
  stringsAsFactors = FALSE
)


# ============================================================
# 7. Keep only Normal and TNBC samples
# ============================================================

sample_map <- sample_map[
  sample_map$Cancer_type %in% c("Normal", "TNBC"),
  ,
  drop = FALSE
]

# Keep only samples present in the count matrix.
sample_map <- sample_map[
  sample_map$ID_REF %in% colnames(counts),
  ,
  drop = FALSE
]


# ============================================================
# 8. Extract candidate-gene expression
# ============================================================

candidate_counts <- counts[
  counts$GeneID %in% candidate_annotation$GeneID,
  ,
  drop = FALSE
]

# Match candidate annotation order.
candidate_counts <- candidate_counts[
  match(
    candidate_annotation$GeneID,
    candidate_counts$GeneID
  ),
  ,
  drop = FALSE
]

# Set gene symbols as row names.
rownames(candidate_counts) <-
  candidate_annotation$Symbol

# Keep only Normal/TNBC sample columns.
sample_ids <- sample_map$ID_REF

candidate_expression <- candidate_counts[
  ,
  sample_ids,
  drop = FALSE
]

# Convert to numeric matrix.
candidate_expression <- as.matrix(
  candidate_expression
)

storage.mode(candidate_expression) <- "numeric"


# ============================================================
# 9. Create discovery ROC results
# ============================================================

discovery_roc_results <- data.frame(
  Gene = character(),
  AUC = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Direction = character(),
  stringsAsFactors = FALSE
)

# Encode TNBC as positive class.
response <- ifelse(
  sample_map$Cancer_type == "TNBC",
  1,
  0
)


# ============================================================
# 10. Calculate ROC/AUC for each candidate gene
# ============================================================

for (gene in rownames(candidate_expression)) {
  
  expression_values <- as.numeric(
    candidate_expression[gene, ]
  )
  
  # Determine direction from the discovery cohort.
  gene_direction <- ifelse(
    mean(
      expression_values[
        sample_map$Cancer_type == "TNBC"
      ]
    ) >
      mean(
        expression_values[
          sample_map$Cancer_type == "Normal"
        ]
      ),
    "Higher in TNBC",
    "Lower in TNBC"
  )
  
  roc_object <- roc(
    response = response,
    predictor = expression_values,
    quiet = TRUE,
    direction = "auto"
  )
  
  auc_value <- as.numeric(
    auc(roc_object)
  )
  
  auc_ci <- as.numeric(
    ci.auc(roc_object)
  )
  
  discovery_roc_results <- rbind(
    discovery_roc_results,
    data.frame(
      Gene = gene,
      AUC = auc_value,
      CI_Lower = auc_ci[1],
      CI_Upper = auc_ci[3],
      Direction = gene_direction,
      stringsAsFactors = FALSE
    )
  )
}


# ============================================================
# 11. Load discovery molecular signatures
# ============================================================

discovery_signatures <- read.csv(
  "Results/Molecular_Signatures/Discovery_Cohort_Signature_Scores.csv",
  stringsAsFactors = FALSE
)

discovery_signatures <- discovery_signatures[
  discovery_signatures$Group %in% c("Normal", "TNBC"),
  ,
  drop = FALSE
]

signature_response <- ifelse(
  discovery_signatures$Group == "TNBC",
  1,
  0
)


# ============================================================
# 12. Calculate signature ROC/AUC
# ============================================================

signature_names <- c(
  "Proliferation_Score",
  "Metabolic_Score"
)

discovery_signature_roc <- data.frame(
  Signature = character(),
  AUC = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Direction = character(),
  stringsAsFactors = FALSE
)

for (signature in signature_names) {
  
  score <- discovery_signatures[[signature]]
  
  roc_object <- roc(
    response = signature_response,
    predictor = score,
    quiet = TRUE,
    direction = "auto"
  )
  
  auc_value <- as.numeric(
    auc(roc_object)
  )
  
  auc_ci <- as.numeric(
    ci.auc(roc_object)
  )
  
  direction <- ifelse(
    mean(
      score[
        discovery_signatures$Group == "TNBC"
      ]
    ) >
      mean(
        score[
          discovery_signatures$Group == "Normal"
        ]
      ),
    "Higher in TNBC",
    "Lower in TNBC"
  )
  
  discovery_signature_roc <- rbind(
    discovery_signature_roc,
    data.frame(
      Signature = signature,
      AUC = auc_value,
      CI_Lower = auc_ci[1],
      CI_Upper = auc_ci[3],
      Direction = direction,
      stringsAsFactors = FALSE
    )
  )
}


# ============================================================
# 13. Load independent molecular signatures
# ============================================================

validation_signatures <- read.csv(
  "Results/Molecular_Signatures/GSE52194_Signature_Scores.csv",
  stringsAsFactors = FALSE
)

validation_signatures <- validation_signatures[
  validation_signatures$Group %in% c("Normal", "TNBC"),
  ,
  drop = FALSE
]

validation_response <- ifelse(
  validation_signatures$Group == "TNBC",
  1,
  0
)


# ============================================================
# 14. Independent validation ROC/AUC
# ============================================================

validation_signature_roc <- data.frame(
  Signature = character(),
  AUC = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Direction = character(),
  stringsAsFactors = FALSE
)

for (signature in signature_names) {
  
  score <- validation_signatures[[signature]]
  
  roc_object <- roc(
    response = validation_response,
    predictor = score,
    quiet = TRUE,
    direction = "auto"
  )
  
  auc_value <- as.numeric(
    auc(roc_object)
  )
  
  auc_ci <- as.numeric(
    ci.auc(roc_object)
  )
  
  direction <- ifelse(
    mean(
      score[
        validation_signatures$Group == "TNBC"
      ]
    ) >
      mean(
        score[
          validation_signatures$Group == "Normal"
        ]
      ),
    "Higher in TNBC",
    "Lower in TNBC"
  )
  
  validation_signature_roc <- rbind(
    validation_signature_roc,
    data.frame(
      Signature = signature,
      AUC = auc_value,
      CI_Lower = auc_ci[1],
      CI_Upper = auc_ci[3],
      Direction = direction,
      stringsAsFactors = FALSE
    )
  )
}

# ============================================================
# 14B. Combine ROC/AUC results
# ============================================================

discovery_gene_summary <- data.frame(
  Analysis = "Individual gene",
  Feature = discovery_roc_results$Gene,
  Cohort = "Discovery",
  AUC = discovery_roc_results$AUC,
  CI_Lower = discovery_roc_results$CI_Lower,
  CI_Upper = discovery_roc_results$CI_Upper,
  Direction = discovery_roc_results$Direction,
  stringsAsFactors = FALSE
)

discovery_signature_summary <- data.frame(
  Analysis = "Molecular signature",
  Feature = discovery_signature_roc$Signature,
  Cohort = "Discovery",
  AUC = discovery_signature_roc$AUC,
  CI_Lower = discovery_signature_roc$CI_Lower,
  CI_Upper = discovery_signature_roc$CI_Upper,
  Direction = discovery_signature_roc$Direction,
  stringsAsFactors = FALSE
)

validation_signature_summary <- data.frame(
  Analysis = "Molecular signature",
  Feature = validation_signature_roc$Signature,
  Cohort = "Independent GSE52194",
  AUC = validation_signature_roc$AUC,
  CI_Lower = validation_signature_roc$CI_Lower,
  CI_Upper = validation_signature_roc$CI_Upper,
  Direction = validation_signature_roc$Direction,
  stringsAsFactors = FALSE
)

roc_auc_summary <- rbind(
  discovery_gene_summary,
  discovery_signature_summary,
  validation_signature_summary
)

# ============================================================
# 15. Create output directories
# ============================================================

if (!dir.exists("Results/ROC_AUC")) {
  
  dir.create(
    "Results/ROC_AUC",
    recursive = TRUE
  )
  
}

if (!dir.exists("Figures/ROC_AUC")) {
  
  dir.create(
    "Figures/ROC_AUC",
    recursive = TRUE
  )
  
}


# ============================================================
# 16. Save ROC/AUC results
# ============================================================

write.csv(
  discovery_roc_results,
  "Results/ROC_AUC/Discovery_Gene_ROC_AUC.csv",
  row.names = FALSE
)

write.csv(
  discovery_signature_roc,
  "Results/ROC_AUC/Discovery_Signature_ROC_AUC.csv",
  row.names = FALSE
)

write.csv(
  validation_signature_roc,
  "Results/ROC_AUC/Independent_Signature_ROC_AUC.csv",
  row.names = FALSE
)

write.csv(
  roc_auc_summary,
  "Results/ROC_AUC/ROC_AUC_Combined_Summary.csv",
  row.names = FALSE
)

# ============================================================
# 17. Create publication-style AUC comparison figure
# ============================================================

plot_data <- data.frame(
  Cohort = c(
    "Discovery",
    "Discovery",
    "Independent GSE52194",
    "Independent GSE52194"
  ),
  Signature = c(
    "Proliferation",
    "Metabolic",
    "Proliferation",
    "Metabolic"
  ),
  AUC = c(
    discovery_signature_roc$AUC[
      discovery_signature_roc$Signature == "Proliferation_Score"
    ],
    discovery_signature_roc$AUC[
      discovery_signature_roc$Signature == "Metabolic_Score"
    ],
    validation_signature_roc$AUC[
      validation_signature_roc$Signature == "Proliferation_Score"
    ],
    validation_signature_roc$AUC[
      validation_signature_roc$Signature == "Metabolic_Score"
    ]
  )
)


# ------------------------------------------------------------
# Plot positions
# ------------------------------------------------------------

x <- c(
  0.72,   # Discovery - Proliferation
  1.72,   # Discovery - Metabolic
  1.28,   # Independent - Proliferation
  2.28    # Independent - Metabolic
)


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

png(
  "Figures/ROC_AUC/Signature_AUC_Comparison.png",
  width = 3200,
  height = 2400,
  res = 300
)


par(
  mar = c(6, 6, 5, 2),
  mgp = c(3.5, 1, 0)
)


# ------------------------------------------------------------
# Main plot
# ------------------------------------------------------------

plot(
  x,
  plot_data$AUC,
  type = "n",
  ylim = c(0.45, 1.10),
  xlim = c(0.45, 2.55),
  xaxt = "n",
  xlab = "Molecular Signature",
  ylab = "Area Under the Curve (AUC)",
  main = "Molecular Signature AUC Comparison",
  cex.main = 1.4,
  cex.lab = 1.2,
  cex.axis = 1.1
)


# ------------------------------------------------------------
# Random-classification reference
# ------------------------------------------------------------

abline(
  h = 0.5,
  lty = 2
)


# ------------------------------------------------------------
# Divider between molecular signatures
# ------------------------------------------------------------

abline(
  v = 1.5,
  lty = 3
)


# ------------------------------------------------------------
# Discovery points
# ------------------------------------------------------------

points(
  x[1:2],
  plot_data$AUC[1:2],
  pch = 19,
  cex = 1.8
)


# ------------------------------------------------------------
# Independent validation points
# ------------------------------------------------------------

points(
  x[3:4],
  plot_data$AUC[3:4],
  pch = 17,
  cex = 1.8
)


# ------------------------------------------------------------
# X-axis labels
# ------------------------------------------------------------

axis(
  1,
  at = c(1, 2),
  labels = c(
    "Proliferation",
    "Metabolic"
  ),
  cex.axis = 1.15
)


# ------------------------------------------------------------
# AUC labels
# ------------------------------------------------------------

# Discovery labels
text(
  x[1],
  1.035,
  "AUC = 1.00",
  cex = 1.0
)

text(
  x[2],
  1.035,
  "AUC = 1.00",
  cex = 1.0
)


# Independent validation labels
text(
  x[3],
  1.035,
  "AUC = 1.00",
  cex = 1.0
)

text(
  x[4],
  1.035,
  "AUC = 1.00",
  cex = 1.0
)


# ------------------------------------------------------------
# Cohort labels
# ------------------------------------------------------------

# Discovery
text(
  x[1],
  0.975,
  "Discovery",
  cex = 1.0
)

text(
  x[2],
  0.975,
  "Discovery",
  cex = 1.0
)


# Independent validation
text(
  x[3],
  0.950,
  "Independent",
  cex = 1.0
)

text(
  x[3],
  0.925,
  "GSE52194",
  cex = 1.0
)

text(
  x[4],
  0.950,
  "Independent",
  cex = 1.0
)

text(
  x[4],
  0.925,
  "GSE52194",
  cex = 1.0
)


# ------------------------------------------------------------
# Footnote
# ------------------------------------------------------------

mtext(
  "Exploratory analysis; independent cohort contains 5 TNBC and 3 Normal samples",
  side = 1,
  line = 4.5,
  cex = 0.9
)


# ------------------------------------------------------------
# Close device
# ------------------------------------------------------------

dev.off()
# ============================================================
# 18. Summary
# ============================================================

cat(
  "ROC/AUC analysis completed.\n"
)

cat(
  "Discovery genes analyzed:",
  nrow(discovery_roc_results),
  "\n"
)

cat(
  "Discovery signatures analyzed:",
  nrow(discovery_signature_roc),
  "\n"
)

cat(
  "Independent signatures analyzed:",
  nrow(validation_signature_roc),
  "\n"
)

cat(
  "Results saved to:",
  "Results/ROC_AUC/",
  "\n"
)

cat(
  "Figure saved to:",
  "Figures/ROC_AUC/Signature_AUC_Comparison.png",
  "\n"
)