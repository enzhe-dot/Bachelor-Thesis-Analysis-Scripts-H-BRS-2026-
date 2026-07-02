# Path to your single count matrix (.tsv)
# Rows = features (miRNA / tRNA), Columns = sample IDs (SRRxxxxxxx)
count_file <- "/on1/MAF/results2/expr_data3/h3.mt_stm_miR.tsv"

# Path to SRA metadata file (tab-separated or comma-separated)
meta_file  <-"/on1/MAF/results2/metadata/SraRunTable_hepato.csv"
meta_sep <- ","

# Output directory (will be created if it does not exist)
OUTPUT_DIR <- "/on1/MAF/results2/output_dir_hepato"

CANCER_NAME = "Hepatocellular"

# Column name in metadata that contains tissue/condition labels
# e.g. "tissue", "source_name", "tissue_type" – check your SraRunTable
CONDITION_COLUMN <- "source_name"

# What string means CANCER in that column?   (partial match, case-insensitive)
CANCER_STRING  <- "tumor"

# What string means NORMAL/CONTROL?          (partial match, case-insensitive)
NORMAL_STRING  <- "normal"

# Column in metadata that holds the Run ID (must match colnames of count matrix)
RUN_COLUMN <- "Run"

# DESeq2 thresholds
LFC_THRESHOLD  <- 1     # |log2FoldChange| cut-off for significance
PADJ_THRESHOLD <- 0.05  # adjusted p-value cut-off

# How many top genes to show in heatmap / stripplot?
TOP_N <- 20


source("/on1/MAF/R_import_scripts/DesrStat_EDA_DESeq2_fixed.R")
