# ============================================================
# CROSS-CANCER ANALYSIS
# Input:  CANCER_RESULTS (list of DESeq2 csv paths)
# ============================================================



# install.packages(c("dplyr", "tibble", "ggplot2", "ggrepel","pheatmap","ggVennDiagram"))



if (!exists("CANCER_RESULTS")) stop("Run via cross_cancer_config.R")

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cancer_names <- names(CANCER_RESULTS)
n_cancers    <- length(cancer_names)
message("Cancers: ", paste(cancer_names, collapse = ", "))

# ── 1. Loading CSV ───────────────────────────────────────────
res_list <- lapply(CANCER_RESULTS, function(path) {
  df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  

  names(df) <- names(df) %>%
    sub("^log2FC$",    "log2FoldChange", .) %>%
    sub("^lfcSE$",     "lfcSE",          .) %>%
    sub("^DE_Status$", "DE_status",      .)
  
  # Убираем строки без LFC или padj
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]
  df
})
names(res_list) <- cancer_names


for (cn in cancer_names) {
  message(cn, " — rows: ", nrow(res_list[[cn]]),
          " | cols: ", paste(colnames(res_list[[cn]]), collapse=", "))
}

# ── 2. Significan genes ──────────────────────────────────────────
sig_list <- lapply(res_list, function(df) {
  df$gene[df$padj < PADJ_THRESHOLD &
            abs(df$log2FoldChange) > LFC_THRESHOLD]
})

sig_up_list <- lapply(res_list, function(df) {
  df$gene[df$padj < PADJ_THRESHOLD &
            df$log2FoldChange > LFC_THRESHOLD]
})

sig_down_list <- lapply(res_list, function(df) {
  df$gene[df$padj < PADJ_THRESHOLD &
            df$log2FoldChange < -LFC_THRESHOLD]
})

for (cn in cancer_names) {
  message(cn, ": ", length(sig_list[[cn]]), " significant (",
          length(sig_up_list[[cn]]),   " up / ",
          length(sig_down_list[[cn]]), " down)")
}


# install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

# ── Venn: all importanat ────────────────────────────────────────
p_venn_all <- ggVennDiagram(sig_list,
                            label_alpha = 0,
                            label       = "count") +
  scale_fill_gradient(low = "white", high = "#E53935") +
  labs(title = "All significant miRNAs") +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "01_venn_all_significant.pdf"),
       p_venn_all, width = 7, height = 6)

# ── Venn: upregulated ─────────────────────────────────────────
p_venn_up <- ggVennDiagram(sig_up_list,
                           label_alpha = 0,
                           label       = "count") +
  scale_fill_gradient(low = "white", high = "#E53935") +
  labs(title = "Upregulated miRNAs") +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "02_venn_upregulated.pdf"),
       p_venn_up, width = 7, height = 6)

# ── Venn: downregulated ───────────────────────────────────────
p_venn_down <- ggVennDiagram(sig_down_list,
                             label_alpha = 0,
                             label       = "count") +
  scale_fill_gradient(low = "white", high = "#1E88E5") +
  labs(title = "Downregulated miRNAs") +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "03_venn_downregulated.pdf"),
       p_venn_down, width = 7, height = 6)

message("  saved -> 01-03 Venn diagrams (ggVennDiagram)")



# ──Specificity table ──────────────────────────────────────
message("\n=== Specificity table ===")

all_genes <- unique(unlist(lapply(res_list, `[[`, "gene")))
spec_mat  <- data.frame(gene = all_genes)

for (cn in cancer_names) {
  spec_mat[[paste0(cn, "_sig")]]  <- spec_mat$gene %in% sig_list[[cn]]
  spec_mat[[paste0(cn, "_dir")]]  <- ifelse(
    spec_mat$gene %in% sig_up_list[[cn]],   "UP",
    ifelse(spec_mat$gene %in% sig_down_list[[cn]], "DOWN", "ns")
  )
 
  lfc_lookup <- setNames(res_list[[cn]]$log2FoldChange,
                         res_list[[cn]]$gene)
  spec_mat[[paste0(cn, "_LFC")]] <- round(
    lfc_lookup[spec_mat$gene], 3)
}

spec_mat$n_cancers_sig <- rowSums(
  spec_mat[, grep("_sig$", colnames(spec_mat))]
)
spec_mat <- spec_mat[spec_mat$n_cancers_sig > 0, ]
spec_mat <- spec_mat[order(-spec_mat$n_cancers_sig), ]

write.table(spec_mat,
            file.path(OUTPUT_DIR, "08_mirna_specificity.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("  saved -> 08_mirna_specificity.tsv")
message("\n=== Cross-cancer analysis DONE ===")
message("Output dir: ", OUTPUT_DIR)
