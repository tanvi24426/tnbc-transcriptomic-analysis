# ============================================================
# 13. TCGA-BRCA Survival Analysis
# ============================================================
#
# Purpose:
# Evaluate the association between expression of the final
# 10 candidate genes and overall survival in TCGA-BRCA.
#
# Survival data:
# GDC TCGA-BRCA
#
# Expression data:
# UCSC Xena GDC TCGA-BRCA STAR-TPM
#
# Analysis:
# - Patient/sample matching
# - Primary tumor sample selection
# - Duplicate primary sample averaging
# - Kaplan-Meier analysis
# - Log-rank test
# - Cox proportional hazards regression
# - Benjamini-Hochberg FDR correction
# - Kaplan-Meier figures
# - Cox forest plot
#
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(jsonlite)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)


# ============================================================
# 2. Candidate genes
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
# 3. Create output directories
# ============================================================

dir.create(
  "Results/Survival",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/Survival",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. Retrieve TCGA-BRCA survival data from GDC
# ============================================================

gdc_survival_url <- paste0(
  "https://api.gdc.cancer.gov/",
  "analysis/survival?",
  "filters=",
  URLencode(
    '{"op":"=","content":{"field":"cases.project.project_id","value":"TCGA-BRCA"}}',
    reserved = TRUE
  ),
  "&pretty=true"
)

gdc_response <- fromJSON(
  gdc_survival_url
)

survival_data <- gdc_response$results$donors[[1]]


# ============================================================
# 5. Validate survival data
# ============================================================

cat(
  "Survival patients:",
  nrow(survival_data),
  "\n"
)

cat(
  "Survival columns:\n"
)

print(
  colnames(survival_data)
)


# ============================================================
# 6. Define survival event
# ============================================================
#
# GDC:
# TRUE  = censored
# FALSE = observed event/death
#
# survival::Surv() requires:
# TRUE  = event
# FALSE = censored
#
# Therefore:
# event = !censored
# ============================================================

survival_data$event <- !survival_data$censored

cat(
  "Deaths:",
  sum(survival_data$event),
  "\n"
)

cat(
  "Censored:",
  sum(!survival_data$event),
  "\n"
)


# ============================================================
# 7. Xena STAR-TPM file
# ============================================================

tpm_file <- file.path(
  "Data",
  "TCGA-BRCA.star_tpm.tsv.gz"
)

if (!file.exists(tpm_file)) {
  
  stop(
    paste(
      "TCGA-BRCA STAR-TPM file not found:",
      tpm_file
    )
  )
  
}


# ============================================================
# 8. GENCODE v36 gene mapping
# ============================================================

probemap_url <- paste0(
  "https://gdc-hub.s3.us-east-1.amazonaws.com/",
  "download/gencode.v36.annotation.gtf.gene.probemap"
)

probemap_file <- tempfile(
  fileext = ".tsv"
)

download.file(
  probemap_url,
  probemap_file,
  mode = "wb"
)

gene_map <- read.delim(
  probemap_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 9. Map candidate genes to Ensembl IDs
# ============================================================

candidate_mapping <- gene_map[
  gene_map$gene %in% candidate_genes,
  ,
  drop = FALSE
]

if (
  nrow(candidate_mapping) !=
  length(candidate_genes)
) {
  
  missing_genes <- setdiff(
    candidate_genes,
    candidate_mapping$gene
  )
  
  stop(
    paste(
      "Candidate genes missing from GENCODE mapping:",
      paste(
        missing_genes,
        collapse = ", "
      )
    )
  )
  
}

candidate_ids <- candidate_mapping$id

cat(
  "Candidate genes mapped:",
  nrow(candidate_mapping),
  "\n"
)


# ============================================================
# 10. Read Xena TPM header
# ============================================================

tpm_header <- read.delim(
  gzfile(tpm_file),
  nrows = 0,
  check.names = FALSE
)

header <- colnames(tpm_header)

cat(
  "Expression samples:",
  length(header) - 1,
  "\n"
)


# ============================================================
# 11. Extract only the 10 candidate genes
# ============================================================
#
# The full matrix contains approximately 60,000 genes.
# We read the compressed file line-by-line so the complete
# matrix does not need to be loaded into memory.
# ============================================================

con <- gzfile(
  tpm_file,
  open = "rt"
)

# Read header
header <- strsplit(
  readLines(con, n = 1),
  "\t",
  fixed = TRUE
)[[1]]

candidate_rows <- list()

while (
  length(
    line <- readLines(
      con,
      n = 1
    )
  ) > 0
) {
  
  fields <- strsplit(
    line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  gene_id <- fields[1]
  
  if (
    gene_id %in%
    candidate_ids
  ) {
    
    candidate_rows[[gene_id]] <- fields
    
  }
  
  if (
    length(
      intersect(
        names(candidate_rows),
        candidate_ids
      )
    ) ==
    length(candidate_ids)
  ) {
    
    break
    
  }
  
}

close(con)


# ============================================================
# 12. Validate extraction
# ============================================================

cat(
  "Candidate genes extracted:",
  length(candidate_rows),
  "\n"
)

if (
  length(candidate_rows) !=
  length(candidate_ids)
) {
  
  missing_ids <- setdiff(
    candidate_ids,
    names(candidate_rows)
  )
  
  stop(
    paste(
      "Candidate genes missing from TPM matrix:",
      paste(
        missing_ids,
        collapse = ", "
      )
    )
  )
  
}


# ============================================================
# 13. Convert extracted rows into expression matrix
# ============================================================

candidate_matrix <- do.call(
  rbind,
  candidate_rows
)

rownames(candidate_matrix) <-
  candidate_matrix[, 1]

candidate_matrix <-
  candidate_matrix[
    ,
    -1,
    drop = FALSE
  ]

candidate_matrix <- apply(
  candidate_matrix,
  2,
  as.numeric
)

rownames(candidate_matrix) <-
  names(candidate_rows)

colnames(candidate_matrix) <-
  header[-1]


# ============================================================
# 14. Map Ensembl IDs back to gene symbols
# ============================================================

gene_symbol_map <- setNames(
  candidate_mapping$gene,
  candidate_mapping$id
)

rownames(candidate_matrix) <-
  gene_symbol_map[
    rownames(candidate_matrix)
  ]


# Validate gene names

if (
  length(
    setdiff(
      candidate_genes,
      rownames(candidate_matrix)
    )
  ) > 0
) {
  
  stop(
    "One or more candidate genes could not be mapped."
  )
  
}

cat(
  "Expression matrix:",
  nrow(candidate_matrix),
  "genes x",
  ncol(candidate_matrix),
  "samples\n"
)


# ============================================================
# 15. Match expression samples to survival patients
# ============================================================

expression_samples <-
  colnames(candidate_matrix)

expression_patient_ids <-
  substr(
    expression_samples,
    1,
    12
  )

expression_patient_map <-
  data.frame(
    sample_id =
      expression_samples,
    patient_id =
      expression_patient_ids,
    stringsAsFactors =
      FALSE
  )

survival_patients <-
  unique(
    survival_data$submitter_id
  )

matched_expression_map <-
  expression_patient_map[
    expression_patient_map$patient_id %in%
      survival_patients,
    ,
    drop = FALSE
  ]

cat(
  "Expression samples:",
  nrow(expression_patient_map),
  "\n"
)

cat(
  "Expression samples with survival data:",
  nrow(matched_expression_map),
  "\n"
)

cat(
  "Unique patients with both datasets:",
  length(
    unique(
      matched_expression_map$patient_id
    )
  ),
  "\n"
)


# ============================================================
# 16. Select primary tumor samples
# ============================================================

primary_expression_map <-
  matched_expression_map[
    grepl(
      "-01",
      matched_expression_map$sample_id
    ),
    ,
    drop = FALSE
  ]

cat(
  "Primary tumor expression samples:",
  nrow(primary_expression_map),
  "\n"
)

cat(
  "Unique primary-tumor patients:",
  length(
    unique(
      primary_expression_map$patient_id
    )
  ),
  "\n"
)


# ============================================================
# 17. Average duplicate primary samples
# ============================================================

primary_samples <-
  primary_expression_map$sample_id

primary_matrix <-
  candidate_matrix[
    ,
    primary_samples,
    drop = FALSE
  ]

primary_patient_ids <-
  primary_expression_map$patient_id


patient_expression <- sapply(
  unique(primary_patient_ids),
  function(patient) {
    
    sample_cols <-
      which(
        primary_patient_ids ==
          patient
      )
    
    if (
      length(sample_cols) == 1
    ) {
      
      primary_matrix[
        ,
        sample_cols
      ]
      
    } else {
      
      rowMeans(
        primary_matrix[
          ,
          sample_cols,
          drop = FALSE
        ],
        na.rm = TRUE
      )
      
    }
    
  }
)

patient_expression <-
  as.matrix(
    patient_expression
  )

rownames(patient_expression) <-
  rownames(primary_matrix)

cat(
  "Patient-level expression matrix:",
  nrow(patient_expression),
  "genes x",
  ncol(patient_expression),
  "patients\n"
)


# ============================================================
# 18. Convert expression matrix to patient-level data frame
# ============================================================

expression_df <-
  as.data.frame(
    t(patient_expression),
    stringsAsFactors = FALSE
  )

expression_df$submitter_id <-
  rownames(
    expression_df
  )


# ============================================================
# 19. Merge expression and survival data
# ============================================================

survival_df <-
  survival_data[
    survival_data$submitter_id %in%
      expression_df$submitter_id,
    c(
      "submitter_id",
      "time",
      "censored",
      "event"
    ),
    drop = FALSE
  ]

survival_merged <-
  merge(
    survival_df,
    expression_df,
    by = "submitter_id"
  )


# ============================================================
# 20. Validate final survival dataset
# ============================================================

cat(
  "Patients in final survival dataset:",
  nrow(survival_merged),
  "\n"
)

cat(
  "Deaths:",
  sum(survival_merged$event),
  "\n"
)

cat(
  "Censored:",
  sum(!survival_merged$event),
  "\n"
)

cat(
  "Missing survival time:",
  sum(
    is.na(
      survival_merged$time
    )
  ),
  "\n"
)

cat(
  "Missing event status:",
  sum(
    is.na(
      survival_merged$event
    )
  ),
  "\n"
)


# ============================================================
# 21. Survival analysis
# ============================================================

survival_results <- list()


for (gene in candidate_genes) {
  
  # --------------------------------------------------------
  # Gene-specific dataset
  # --------------------------------------------------------
  
  gene_data <-
    survival_merged[
      ,
      c(
        "submitter_id",
        "time",
        "event",
        gene
      )
    ]
  
  colnames(gene_data)[4] <-
    "expression"
  
  gene_data <-
    gene_data[
      complete.cases(
        gene_data
      ),
      ,
      drop = FALSE
    ]
  
  
  # --------------------------------------------------------
  # Median expression split
  # --------------------------------------------------------
  
  median_expression <-
    median(
      gene_data$expression,
      na.rm = TRUE
    )
  
  gene_data$expression_group <-
    ifelse(
      gene_data$expression >=
        median_expression,
      "High",
      "Low"
    )
  
  gene_data$expression_group <-
    factor(
      gene_data$expression_group,
      levels = c(
        "Low",
        "High"
      )
    )
  
  
  # --------------------------------------------------------
  # Kaplan-Meier model
  # --------------------------------------------------------
  
  km_fit <-
    survfit(
      Surv(
        time,
        event
      ) ~ expression_group,
      data = gene_data
    )
  
  
  # --------------------------------------------------------
  # Log-rank test
  # --------------------------------------------------------
  
  logrank_test <-
    survdiff(
      Surv(
        time,
        event
      ) ~ expression_group,
      data = gene_data
    )
  
  logrank_p <-
    1 -
    pchisq(
      logrank_test$chisq,
      df = 1
    )
  
  
  # --------------------------------------------------------
  # Cox proportional hazards model
  # --------------------------------------------------------
  
  cox_fit <-
    coxph(
      Surv(
        time,
        event
      ) ~ expression_group,
      data = gene_data
    )
  
  cox_summary <-
    summary(
      cox_fit
    )
  
  hazard_ratio <-
    cox_summary$coefficients[
      "expression_groupHigh",
      "exp(coef)"
    ]
  
  ci_lower <-
    cox_summary$conf.int[
      "expression_groupHigh",
      "lower .95"
    ]
  
  ci_upper <-
    cox_summary$conf.int[
      "expression_groupHigh",
      "upper .95"
    ]
  
  cox_p <-
    cox_summary$coefficients[
      "expression_groupHigh",
      "Pr(>|z|)"
    ]
  
  
  # --------------------------------------------------------
  # Store results
  # --------------------------------------------------------
  
  survival_results[[gene]] <-
    data.frame(
      Gene = gene,
      N = nrow(gene_data),
      High_N =
        sum(
          gene_data$expression_group ==
            "High"
        ),
      Low_N =
        sum(
          gene_data$expression_group ==
            "Low"
        ),
      Deaths =
        sum(
          gene_data$event
        ),
      Median_Expression =
        median_expression,
      LogRank_P =
        logrank_p,
      Hazard_Ratio =
        hazard_ratio,
      CI_Lower =
        ci_lower,
      CI_Upper =
        ci_upper,
      Cox_P =
        cox_p,
      stringsAsFactors =
        FALSE
    )
  
  
  # --------------------------------------------------------
  # Kaplan-Meier figure
  # --------------------------------------------------------
  
  km_plot <-
    ggsurvplot(
      km_fit,
      data = gene_data,
      pval = TRUE,
      risk.table = TRUE,
      conf.int = FALSE,
      xlab = "Time (days)",
      ylab =
        "Overall Survival Probability",
      title =
        paste(
          gene,
          "Expression and Overall Survival"
        ),
      legend.title =
        "Expression",
      legend.labs =
        c(
          "Low",
          "High"
        ),
      risk.table.height =
        0.25,
      ggtheme =
        theme_minimal()
    )
  
  
  # --------------------------------------------------------
  # Save KM plot
  # --------------------------------------------------------
  
  ggsave(
    filename =
      paste0(
        "Figures/Survival/",
        gene,
        "_Kaplan_Meier.png"
      ),
    plot =
      km_plot$plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # --------------------------------------------------------
  # Save KM plot with risk table
  # --------------------------------------------------------
  
  png(
    filename =
      paste0(
        "Figures/Survival/",
        gene,
        "_Kaplan_Meier_RiskTable.png"
      ),
    width = 2400,
    height = 2400,
    res = 300
  )
  
  print(km_plot)
  
  dev.off()
  
}


# ============================================================
# 22. Combine survival results
# ============================================================

survival_results_table <-
  bind_rows(
    survival_results
  )


# ============================================================
# 23. Multiple-testing correction
# ============================================================

survival_results_table$LogRank_FDR <-
  p.adjust(
    survival_results_table$LogRank_P,
    method = "BH"
  )

survival_results_table$Cox_FDR <-
  p.adjust(
    survival_results_table$Cox_P,
    method = "BH"
  )


# ============================================================
# 24. Save survival results
# ============================================================

write.csv(
  survival_results_table,
  "Results/Survival/10_Gene_Survival_Analysis.csv",
  row.names = FALSE
)


# ============================================================
# 25. Final Cox forest plot
# ============================================================

forest_data <-
  survival_results_table

forest_data$Gene <-
  factor(
    forest_data$Gene,
    levels =
      forest_data$Gene[
        order(
          forest_data$Hazard_Ratio
        )
      ]
  )

forest_data$Highlight <-
  ifelse(
    forest_data$Cox_P < 0.05,
    "Nominally significant",
    "Not significant"
  )


forest_plot <-
  ggplot(
    forest_data,
    aes(
      x = Hazard_Ratio,
      y = Gene
    )
  ) +
  
  # HR = 1 reference
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey45"
  ) +
  
  # 95% CI
  geom_errorbar(
    aes(
      xmin = CI_Lower,
      xmax = CI_Upper,
      color = Highlight
    ),
    orientation = "y",
    height = 0.18,
    linewidth = 0.9
  ) +
  
  # HR point
  geom_point(
    aes(
      color = Highlight
    ),
    size = 4
  ) +
  
  # HR label
  geom_text(
    aes(
      label =
        sprintf(
          "HR %.2f",
          Hazard_Ratio
        )
    ),
    nudge_y = 0.22,
    size = 3.5,
    color = "black"
  ) +
  
  scale_color_manual(
    values = c(
      "Nominally significant" =
        "#C62828",
      "Not significant" =
        "grey35"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    name =
      "Hazard Ratio (95% CI)",
    breaks = c(
      0.5,
      0.75,
      1,
      1.25,
      1.5,
      1.75,
      2
    ),
    limits = c(
      0.45,
      2.15
    )
  ) +
  
  labs(
    title =
      "Overall Survival Association of 10 Candidate Genes",
    subtitle =
      paste0(
        "TCGA-BRCA cohort (n = ",
        unique(
          forest_data$N
        ),
        "; deaths = ",
        unique(
          forest_data$Deaths
        ),
        ")"
      ),
    y = NULL,
    caption =
      paste0(
        "Hazard ratios >1 indicate higher mortality risk with high gene expression.\n",
        "PLK1: Cox p = 0.021; FDR = 0.214. ",
        "Nominal association did not remain significant after multiple-testing correction."
      )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 17
      ),
    plot.subtitle =
      element_text(
        size = 11
      ),
    axis.title.x =
      element_text(
        face = "bold",
        size = 13
      ),
    axis.text.y =
      element_text(
        face = "bold",
        size = 12
      ),
    axis.text.x =
      element_text(
        size = 11
      ),
    legend.position =
      "none",
    plot.caption =
      element_text(
        size = 9.5,
        hjust = 0
      ),
    panel.grid.major.x =
      element_line(
        color = "grey90"
      ),
    panel.grid.major.y =
      element_blank(),
    plot.margin =
      margin(
        15,
        20,
        15,
        15
      )
  )


# Display forest plot
print(
  forest_plot
)


# Save final forest plot
ggsave(
  "Figures/Survival/10_Gene_Cox_Forest_Plot_Final.png",
  forest_plot,
  width = 9,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 26. Final summary
# ============================================================

cat(
  "\n",
  "============================================================\n",
  "Survival analysis completed.\n",
  "============================================================\n",
  "\n"
)

cat(
  "Candidate genes:",
  length(candidate_genes),
  "\n"
)

cat(
  "Final patients:",
  nrow(survival_merged),
  "\n"
)

cat(
  "Deaths:",
  sum(survival_merged$event),
  "\n"
)

cat(
  "Censored:",
  sum(!survival_merged$event),
  "\n"
)

cat(
  "Results saved to:",
  "Results/Survival/\n"
)

cat(
  "Forest plot saved to:",
  "Figures/Survival/10_Gene_Cox_Forest_Plot_Final.png\n"
)

cat(
  "\nTop nominal survival association:\n"
)

print(
  survival_results_table[
    order(
      survival_results_table$Cox_P
    ),
    c(
      "Gene",
      "Hazard_Ratio",
      "CI_Lower",
      "CI_Upper",
      "Cox_P",
      "Cox_FDR"
    )
  ]
)

cat(
  "\n============================================================\n"
)