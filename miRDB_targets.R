library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ---- 1. Paths ----
mirdb_gz <- "/on1/MAF/data/miRDB_v6.0_prediction_result.txt"

mirna_lists <- c(
  Breast         = "/on1/MAF/results2/breast_significant_miRNA.csv",
  Lung           = "/on1/MAF/results2/lung_significan_miRNA.csv",
  Larynx         = "/on1/MAF/results2/laryngeal_significant_miRNA.csv",
  Hepatocellular = "/on1/MAF/results2/hepatocellular_significant_miR.csv"
)

output_dir <- "/on1/MAF/results2/targets MirDB"
dir.create(output_dir, showWarnings = FALSE)


# ---- 2. Reading miRDB ----
cat("Читаю miRDB v6.0...\n")
con <- gzfile(mirdb_gz, "rt")
mirdb <- read.table(con, sep = "\t", header = FALSE,
                    col.names = c("miRNA", "RefSeq", "Score"),
                    quote = "", comment.char = "",
                    stringsAsFactors = FALSE)
close(con)

cat(sprintf("Overall listings: %d\n", nrow(mirdb)))

mirdb_hsa <- mirdb %>%
  filter(grepl("^hsa-", miRNA), Score >= 80)

cat(sprintf("Human, score >= 80: %d\n", nrow(mirdb_hsa)))


# ---- 3. Helper function to read a single miRNA file ----
# Returns a character vector in the original order, without the header.
#
# FIX #1: Previously, the pattern ^(gene|miRNA|mirna|miR|V1) could accidentally
# match a real miRNA like "hsa-miR-21-5p" (due to the "miR" substring).
# Now, we remove the line as a header ONLY if it does NOT start
# with "hsa-" (i.e., it is clearly not a genuine miRNA name).

read_mirna_list <- function(path, cancer_name) {
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  
  # for (Excel)
  raw <- gsub("^\ufeff", "", raw)
  
 
  raw <- trimws(raw)
  raw <- raw[nchar(raw) > 0]
  
  # CSV 
  if (any(grepl(",", raw))) {
    raw <- sub(",.*$", "", raw)
    raw <- trimws(raw)
  }
  
  # Fix
  first <- raw[1]
  is_real_mirna <- grepl("^[a-z]{3}-", first, ignore.case = TRUE)  # hsa-, mmu-, ...
  if (!is_real_mirna) {
    cat(sprintf("  [%s] removing title: '%s'\n", cancer_name, first))
    raw <- raw[-1]
  }
  
  raw <- raw[nchar(raw) > 0]
  cat(sprintf("  [%s] %d miRNA, first: '%s'\n", cancer_name, length(raw), raw[1]))
  return(raw)
}


# ---- 4. Reading all lists and preserving the order of miRNA----


cancer_mirna_lists <- lapply(names(mirna_lists), function(cn) {
  read_mirna_list(mirna_lists[cn], cn)
})
names(cancer_mirna_lists) <- names(mirna_lists)


# ---- 5. Mapping RefSeq → Gene Symbol ----
all_mirnas_unique <- unique(unlist(cancer_mirna_lists))
cat(sprintf("\nOverall unique miRNA among all cancers: %d\n", length(all_mirnas_unique)))


mirdb_filtered <- mirdb_hsa %>%
  filter(miRNA %in% all_mirnas_unique)

cat(sprintf("miRNA-target before annotation: %d\n", nrow(mirdb_filtered)))


mirdb_filtered$RefSeq_clean <- gsub("\\..*$", "", mirdb_filtered$RefSeq)

# mapping
refseq2eg <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(mirdb_filtered$RefSeq_clean),
  columns = c("ENTREZID", "SYMBOL", "GENENAME"),
  keytype = "REFSEQ"
)

# adding annotation
mirdb_annot <- mirdb_filtered %>%
  left_join(refseq2eg, by = c("RefSeq_clean" = "REFSEQ")) %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::select(miRNA, Gene_Symbol = SYMBOL, Score,
                Gene_Description = GENENAME, Entrez_ID = ENTREZID,
                RefSeq)

cat(sprintf("miRNA-target пар после аннотации: %d\n", nrow(mirdb_annot)))


# ---- 6. Function building a table for each miRNA----


build_cancer_table <- function(mirna_vec, cancer_name, annot_df) {
  
  order_df <- data.frame(
    miRNA     = mirna_vec,
    row_order = seq_along(mirna_vec),
    stringsAsFactors = FALSE
  )
  # Order cheek 
  result <- order_df %>%
    left_join(annot_df, by = "miRNA") %>%
    filter(!is.na(Gene_Symbol)) %>%
   
    arrange(row_order, desc(Score)) %>%
    dplyr::select(miRNA, Gene_Symbol, Score, Gene_Description, Entrez_ID, RefSeq)
  
  # which mRNA were not found 
  not_found <- setdiff(mirna_vec, unique(result$miRNA))
  if (length(not_found) > 0) {
    cat(sprintf("  [%s] not fpund in  miRDB (score >= 80): %s\n",
                cancer_name, paste(not_found, collapse = ", ")))
  }
  
  cat(sprintf("  [%s] %d miRNA-target пар, %d уникальных miRNA\n",
              cancer_name, nrow(result), length(unique(result$miRNA))))
  return(result)
}


# ---- 7.Build and safe ----
cat("\nI am making a table...\n")

all_cancer_tables <- list()

for (cn in names(cancer_mirna_lists)) {
  tbl <- build_cancer_table(cancer_mirna_lists[[cn]], cn, mirdb_annot)
  all_cancer_tables[[cn]] <- tbl
  
  out_path <- file.path(output_dir, paste0("miRDB_targets_", cn, ".csv"))
  write.csv(tbl, out_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", out_path))
}

#combined table 
combined <- bind_rows(
  lapply(names(all_cancer_tables), function(cn) {
    all_cancer_tables[[cn]] %>% mutate(Cancer = cn)
  })
) %>% dplyr::select(Cancer, everything())

out_all <- file.path(output_dir, "miRDB_targets_all.csv")
write.csv(combined, out_all, row.names = FALSE)
cat(sprintf("\nCommon table: %s\n", out_all))


# ---- 8. Outcome ----
cat("\n===== Outcome =====\n")
for (cn in names(all_cancer_tables)) {
  tbl <- all_cancer_tables[[cn]]
  cat(sprintf("[%s] miRNA: %d | Gene targets: %d | пар: %d\n",
              cn,
              length(unique(tbl$miRNA)),
              length(unique(tbl$Gene_Symbol)),
              nrow(tbl)))
}

cat("\nFirst 20 lines  (Breast):\n")
print(head(all_cancer_tables[["Breast"]], 20))

