# ============================================================
# 01_explore_data.R
# Data Exploration and Sample Preparation
# ============================================================


# ------------------------------------------------------------
# 1. Check project structure
# ------------------------------------------------------------

# Display the current working directory
getwd()

# Display files in the Data folder
list.files("Data")


# ------------------------------------------------------------
# 2. Load GEO sample metadata
# ------------------------------------------------------------

# Load the mapping between GEO sample IDs and
# biological sample information
sample_map <- read.csv(
  "Data/GSE233242_Sample_ID_to_GEO_ids.csv.gz",
  stringsAsFactors = FALSE
)

# Inspect the metadata
head(sample_map)

# Check the total number of samples
nrow(sample_map)

# Display the available cancer types
table(sample_map$Cancer_type)


# ------------------------------------------------------------
# 3. Identify Normal and TNBC samples
# ------------------------------------------------------------

# Select Normal and TNBC samples for the analysis
analysis_metadata <- sample_map[
  sample_map$Cancer_type %in% c(
    "Normal",
    "TNBC"
  ),
]

# Check the number of samples in each group
table(analysis_metadata$Cancer_type)


# ------------------------------------------------------------
# 4. Identify matched patients
# ------------------------------------------------------------

# Identify patients for whom both Normal and TNBC
# samples are available
paired_patients <- names(
  which(
    table(analysis_metadata$Patient_ID) >= 2
  )
)

# Display the paired patient IDs
paired_patients

# Restrict the metadata to matched patients
analysis_metadata <- analysis_metadata[
  analysis_metadata$Patient_ID %in% paired_patients,
]

# Check the sample distribution after patient matching
table(analysis_metadata$Cancer_type)


# ------------------------------------------------------------
# 5. Load raw RNA-seq count data
# ------------------------------------------------------------

# Load the raw gene-level count matrix
counts <- read.delim(
  "Data/GSE233242_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  check.names = FALSE
)

# Inspect the dimensions of the count matrix
dim(counts)

# Inspect the first few rows
head(counts)


# ------------------------------------------------------------
# 6. Check sample ID matching
# ------------------------------------------------------------

# Extract sample IDs from the metadata
metadata_sample_ids <- analysis_metadata$ID_REF

# Identify metadata samples that are present
# in the raw count matrix
matched_sample_ids <- metadata_sample_ids[
  metadata_sample_ids %in% colnames(counts)
]

# Display samples that are present in the metadata
# but missing from the count matrix
metadata_sample_ids[
  !metadata_sample_ids %in% colnames(counts)
]

# Count successfully matched samples
length(matched_sample_ids)


# ------------------------------------------------------------
# 7. Define the final analysis sample set
# ------------------------------------------------------------

# The original analysis used eight matched
# Normal-TNBC patient pairs, giving 16 samples.
#
# GSM7416157 is absent from the raw count matrix.
# Therefore, its paired TNBC sample GSM7416158
# is excluded to maintain complete patient pairs.

analysis_sample_ids <- c(
  "GSM7416122",
  "GSM7416123",
  "GSM7416165",
  "GSM7416166",
  "GSM7416154",
  "GSM7416156",
  "GSM7416139",
  "GSM7416140",
  "GSM7416137",
  "GSM7416138",
  "GSM7416172",
  "GSM7416173",
  "GSM7416174",
  "GSM7416175",
  "GSM7416111",
  "GSM7416112"
)

# Verify that all selected samples are present
# in the raw count matrix
analysis_sample_ids %in% colnames(counts)

# Count the successfully matched analysis samples
sum(
  analysis_sample_ids %in% colnames(counts)
)


# ------------------------------------------------------------
# 8. Create the final analysis metadata
# ------------------------------------------------------------

# Keep only the selected 16 samples
analysis_metadata <- analysis_metadata[
  analysis_metadata$ID_REF %in% analysis_sample_ids,
]

# Reorder metadata to match the order of
# analysis_sample_ids
analysis_metadata <- analysis_metadata[
  match(
    analysis_sample_ids,
    analysis_metadata$ID_REF
  ),
]

# Confirm the final sample distribution
table(
  analysis_metadata$Cancer_type
)

# Confirm the final number of samples
nrow(analysis_metadata)


# ------------------------------------------------------------
# 9. Create the final count matrix
# ------------------------------------------------------------

# Keep GeneID and the 16 selected analysis samples
counts_subset <- counts[
  ,
  c(
    "GeneID",
    analysis_sample_ids
  )
]

# Confirm the dimensions of the selected matrix
dim(counts_subset)

# Display the selected sample IDs
colnames(counts_subset)[-1]


# ------------------------------------------------------------
# 10. Inspect gene identifiers
# ------------------------------------------------------------

# Check the first few GeneIDs
head(
  counts_subset$GeneID
)

# Check whether GeneIDs are duplicated
sum(
  duplicated(counts_subset$GeneID)
)


# ------------------------------------------------------------
# 11. Calculate sequencing depth
# ------------------------------------------------------------

# Calculate total read counts for each sample
library_sizes <- colSums(
  counts_subset[, -1]
)

# Display sequencing depth
library_sizes

# Create a basic sequencing-depth plot
barplot(
  library_sizes,
  las = 2,
  main = "Sequencing Depth Across Samples",
  ylab = "Total Counts"
)


# ------------------------------------------------------------
# 12. Examine count distributions
# ------------------------------------------------------------

# Calculate log2-transformed counts
# for visualization of count distributions
log_counts <- log2(
  counts_subset[, -1] + 1
)

# Plot the distribution of log-transformed counts
boxplot(
  log_counts,
  las = 2,
  main = "Log2 Count Distribution",
  ylab = "log2(count + 1)"
)


# ------------------------------------------------------------
# 13. Save basic QC figures
# ------------------------------------------------------------

# Create the Figures folder if it does not exist
if (!dir.exists("Figures")) {
  dir.create("Figures")
}


# Save sequencing-depth plot
png(
  "Figures/sequencing_depth.png",
  width = 1800,
  height = 1400,
  res = 250
)

barplot(
  library_sizes,
  las = 2,
  main = "Sequencing Depth Across Samples",
  ylab = "Total Counts"
)

dev.off()


# Save count-distribution plot
png(
  "Figures/gene_count_distribution.png",
  width = 1800,
  height = 1400,
  res = 250
)

boxplot(
  log_counts,
  las = 2,
  main = "Log2 Count Distribution",
  ylab = "log2(count + 1)"
)

dev.off()


# ------------------------------------------------------------
# 14. Final checks
# ------------------------------------------------------------

# Confirm the final number of samples
length(analysis_sample_ids)

# Confirm the Normal/TNBC distribution
table(
  analysis_metadata$Cancer_type
)

# Confirm the count matrix dimensions
dim(counts_subset)

# Confirm that all selected samples are present
all(
  analysis_sample_ids %in% colnames(counts)
)