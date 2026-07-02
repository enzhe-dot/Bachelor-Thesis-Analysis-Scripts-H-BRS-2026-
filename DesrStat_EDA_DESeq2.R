
#    BiocManager::install(c("DESeq2", "apeglm"))
# install.packages(c("pheatmap", "RColorBrewer", "ggplot2", "ggrepel","gghalves"))



suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

  # loads: count_file, meta_file, OUTPUT_DIR,
#        CONDITION_COLUMN, CANCER_STRING, NORMAL_STRING,
#        RUN_COLUMN, LFC_THRESHOLD, PADJ_THRESHOLD, TOP_N

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# helper: save a ggplot to PDF
save_gg <- function(plot_obj, filename, w = 8, h = 6, limitsize = TRUE) {
  out <- file.path(OUTPUT_DIR, paste0(filename, ".pdf"))
  ggsave(out, plot_obj, width = w, height = h, limitsize = limitsize)
  message("  saved -> ", filename, ".pdf")
}


# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
message("\n=== 1. Loading data ===")

# --- Count matrix ---
# Expected format: first column = feature names, remaining columns = samples (SRRxxxxxx)
count_data <- read.delim(count_file, header = TRUE,
                         check.names = FALSE, sep = "\t")


rownames(count_data) <- count_data[, 1]
count_data           <- count_data[, -1]          # removing first column
count_data           <- round(as.matrix(count_data))
storage.mode(count_data) <- "integer"

message("  Count matrix: ", nrow(count_data), " features x ", ncol(count_data), " samples")

# --- Metadata ---
# Try tab first, then comma
meta_data <- tryCatch(
  read.delim(meta_file, header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, sep = meta_sep),
  error = function(e)
    read.csv(meta_file, header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE)
)
ncol(meta_data)
head(meta_data$source_name)
# loock up
head(colnames(count_data), 3)
head(meta_data$Run, 3)


chartr("", "", colnames(count_data)[1])   
nchar(colnames(count_data)[1])
nchar(meta_data$Run[1])
unique(meta_data$source_name)
nchar(unique(meta_data$source_name))

message("  Metadata: ", nrow(meta_data), " rows, ",
        ncol(meta_data), " columns",nchar(unique(meta_data$source_name)) )

# check required columns exist
if (!CONDITION_COLUMN %in% names(meta_data))
  stop("Column '", CONDITION_COLUMN, "' not found in metadata.\n",
       "Available columns: ", paste(names(meta_data), collapse = ", "))
if (!RUN_COLUMN %in% names(meta_data))
  stop("Column '", RUN_COLUMN, "' not found in metadata.")


# ── 2. PREPARE METADATA ───────────────────────────────────────────────────────
message("\n=== 2. Preparing metadata ===")

meta_data$condition <- ifelse(
  grepl(CANCER_STRING, meta_data[[CONDITION_COLUMN]], ignore.case = TRUE),
  "Cancer",
  ifelse(
    grepl(NORMAL_STRING, meta_data[[CONDITION_COLUMN]], ignore.case = TRUE),
    "Normal",
    NA
  )
)

# Drop samples with unrecognised condition
n_before <- nrow(meta_data)
meta_data <- meta_data[!is.na(meta_data$condition), ]
message("  Dropped ", n_before - nrow(meta_data),
        " rows with unrecognised condition label")

meta_data$condition <- factor(meta_data$condition,
                              levels = c("Normal", "Cancer"))

rownames(meta_data) <- meta_data[[RUN_COLUMN]]

message("  Cancer samples: ", sum(meta_data$condition == "Cancer"),
        "   Normal samples: ", sum(meta_data$condition == "Normal"))


# ── 3. ALIGN COUNT MATRIX & METADATA ─────────────────────────────────────────
message("\n=== 3. Aligning samples ===")

common_samples <- intersect(colnames(count_data), rownames(meta_data))

if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between count matrix and metadata!\n",
       "Count matrix columns (first 5): ",
       paste(head(colnames(count_data), 5), collapse = ", "), "\n",
       "Metadata Run IDs (first 5): ",
       paste(head(rownames(meta_data), 5), collapse = ", "))
}

count_data_filtered <- count_data[, common_samples, drop = FALSE]
meta_data_filtered  <- meta_data[common_samples, , drop = FALSE]


meta_data_filtered <- meta_data_filtered[colnames(count_data_filtered), ]
stopifnot(all(colnames(count_data_filtered) == rownames(meta_data_filtered)))

message("  Common samples used: ", length(common_samples))


# ── 4. EDA  –  RAW DATA QC ────────────────────────────────────────────────────
message("\n=== 4. EDA – Raw QC ===")

## 4.1  Library sizes -----------------------------------------------------------
lib_sizes <- colSums(count_data_filtered)
lib_df    <- data.frame(
  Sample    = factor(names(lib_sizes),
                     levels = names(sort(lib_sizes))),   # sorted order
  TotalReads = lib_sizes,
  Condition  = meta_data_filtered$condition
)

p_lib <- ggplot(lib_df, aes(x = Sample, y = TotalReads / 1e6,
                            fill = Condition)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(Normal = "#43A047", Cancer = "#E53935")) +
  labs(title = "Library Size per Sample",
       x = NULL, y = "Total Reads (Millions)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 7))
save_gg(p_lib, "01_library_sizes",
        h = max(4, length(lib_sizes) * 0.28))





## 4.2  Print summary table
sparsity_pct <- round(100 * mean(count_data_filtered == 0), 1)
cat("\n--- Descriptive Summary ---\n")
cat("Features          :", nrow(count_data_filtered), "\n")
cat("Samples           :", ncol(count_data_filtered), "\n")
cat("Cancer samples    :", sum(meta_data_filtered$condition == "Cancer"), "\n")
cat("Normal samples    :", sum(meta_data_filtered$condition == "Normal"), "\n")
cat("Median lib size   :", round(median(lib_sizes) / 1e6, 2), "M reads\n")
cat("Median detected   :", round(median(det_df$Detected)), "features\n")
cat("Sparsity (% zeros):", sparsity_pct, "%\n\n")


summary_df <- data.frame(
  Metric = c(
    "Total features",
    "Total samples",
    "Cancer samples",
    "Normal samples",
    "Median library size (M reads)",
    "Median detected features per sample",
    "Sparsity (% zeros in matrix)"
  ),
  Value = c(
    nrow(count_data_filtered),
    ncol(count_data_filtered),
    sum(meta_data_filtered$condition == "Cancer"),
    sum(meta_data_filtered$condition == "Normal"),
    round(median(lib_sizes) / 1e6, 2),
    round(median(det_df$Detected)),
    sparsity_pct
  )
)

write.csv(summary_df,
          file = file.path(OUTPUT_DIR, "00_descriptive_summary.csv"),
          row.names = FALSE)

message("  saved -> 00_descriptive_summary.csv")


# ── 5. BUILD DESeqDataSet ─────────────────────────────────────────────────────
message("=== 5. Building DESeqDataSet ===")

meta_clean <- data.frame(
  condition = meta_data_filtered$condition,
  row.names = rownames(meta_data_filtered)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData   = meta_clean,
  design    = ~ condition
)

# Pre-filter: keep features with >= 10 reads in at least 20% of samples
# (slightly lenient for miRNA which can have very focused expression)


keep_samples <- colSums(counts(dds)) > 0
dds <- dds[, keep_samples]
message("Removed zero samples: ", sum(!keep_samples))
# left for breast 7 - normal; 53 - cancer


keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]
message("Genes left : ", sum(keep_genes))


dds <- estimateSizeFactors(dds, type = "poscounts")

message("  Features after filtering: ", nrow(dds),
        "  (removed ", sum(!keep_genes), ")")



# ── 6. NORMALISATION QC ───────────────────────────────────────────────────────
message("\n=== 6. Normalisation QC ===")

## VST for visualisation (blind = TRUE = unbiased by design)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- assay(vsd)
colnames(vst_matrix) <- paste0(CANCER_NAME, "_", colnames(vst_matrix))

write.csv(vst_matrix,
          file = file.path(OUTPUT_DIR, paste0(CANCER_NAME, "_vst.csv")))

meta_export <- data.frame(
  sample    = colnames(vst_matrix),
  cancer    = CANCER_NAME,
  condition = dds$condition,
  row.names = NULL
)
write.csv(meta_export,
          file = file.path(OUTPUT_DIR, paste0(CANCER_NAME, "_meta.csv")),
          row.names = FALSE)

## 6.1  Mean vs. Variance (raw normalised counts, log10 axes) ------------------
# Run DESeq2 now so we have size-factor-normalised counts for the MV plot
# (DESeq() also estimates dispersions needed for plotDispEsts)
dds <- DESeq(dds)

gene_means <- rowMeans(counts(dds, normalized = TRUE))
gene_vars  <- apply(counts(dds, normalized = TRUE), 1, var)   

mv_df <- data.frame(Mean = gene_means, Var = gene_vars)
mv_df <- mv_df[mv_df$Mean > 0 & mv_df$Var > 0, ]

p_mv <- ggplot(mv_df, aes(x = Mean, y = Var)) +
  geom_point(alpha = 0.25, size = 0.7, colour = "#37474F") +
  geom_abline(slope = 1, intercept = 0,
              colour = "red", linetype = "dashed", linewidth = 0.9) +
  scale_x_log10() + scale_y_log10() +
  labs(title    = "Mean-Variance Relationship (normalised counts)",
       subtitle = "Red dashed = Poisson (var = mean)\nPoints above line = overdispersion -> DESeq2 NB model is correct",
       x = "Mean (log10)", y = "Variance (log10)") +
  theme_bw(base_size = 12)
save_gg(p_mv, "04_mean_variance")



# ── Cook's distance plot ──────────────────────────────────────────────────────
message("  Cook's distance plot")

cooks_mat <- assays(dds)[["cooks"]]

# Boxplot Cook's distance по сэмплам
pdf(file.path(OUTPUT_DIR, "C4_cooks_distance.pdf"), width = 12, height = 5)
par(mar = c(8, 5, 2, 2))
boxplot(
  log10(cooks_mat),
  las     = 2,                          
  cex.axis = 0.6,
  main    = "Cook's Distance per Sample (log10)",
  ylab    = "log10(Cook's distance)",
  col     = ifelse(colnames(cooks_mat) %in%
                     colnames(dds)[dds$condition == "Cancer"],
                   "#E53935", "#43A047")
)
abline(h = log10(0.01), col = "blue",  lty = 2)   # threshold 
abline(h = log10(0.1),  col = "red",   lty = 2)
legend("bottomright",
       legend = c("threshold 0.01", "threshold 0.1", "Cancer", "Normal"),
       col    = c("blue", "red", "#E53935", "#43A047"),
       lty    = c(2, 2, NA, NA),
       pch    = c(NA, NA, 15, 15),
       cex = 0.5,
       )

dev.off()
message("  saved -> C4_cooks_distance.pdf")

# suspicous high Cook's distance
max_cooks_per_sample <- apply(cooks_mat, 2, max, na.rm = TRUE)
suspicious_samples   <- names(max_cooks_per_sample)[
  max_cooks_per_sample > quantile(max_cooks_per_sample, 0.95)
]
cat("\nSamples with high Cook's distance (higher that other 95% of samples):\n")
print(round(sort(max_cooks_per_sample[suspicious_samples], 
                 decreasing = TRUE), 3))


## 6.2  Dispersion estimates ----------------------------------------------------
pdf(file.path(OUTPUT_DIR, "05_dispersion_estimates.pdf"), width = 7, height = 6)
plotDispEsts(dds,
             main = paste0("Dispersion Estimates\n",
                           "black = per-gene  |  red = fitted  |  blue = shrunken"))
dev.off()
message("  saved -> 05_dispersion_estimates.pdf")


# ── 7. EDA  –  SAMPLE RELATIONSHIPS ──────────────────────────────────────────
message("\n=== 7. EDA – Sample relationships ===")

## 7.1  PCA --------------------------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data,
                aes(x = PC1, y = PC2, colour = condition, label = name)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 25) +
  scale_colour_manual(values = c(Normal = "#43A047", Cancer = "#E53935")) +
  labs(title    = "PCA – VST-normalised counts",
       subtitle = paste0("Each dot = one sample  |  ",
                         "Separation = different expression profiles"),
       x = paste0("PC1  (", pct_var[1], "% variance)"),
       y = paste0("PC2  (", pct_var[2], "% variance)"),
       colour = "Condition") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")
save_gg(p_pca, "06_PCA", w = 7, h = 6)

## 7.2  Sample-distance heatmap ------------------------------------------------
sampleDists   <- dist(t(assay(vsd)))
sampleDistMat <- as.matrix(sampleDists)
colors_hm     <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

ann_col <- data.frame(Condition = meta_data_filtered$condition,
                      row.names = rownames(meta_data_filtered))

pdf(file.path(OUTPUT_DIR, "07_sample_distance_heatmap.pdf"), width = 9, height = 8)
pheatmap(sampleDistMat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = ann_col,
         col  = colors_hm,
         main = "Sample-to-Sample Distance (Euclidean, VST)\nOutliers = lone bright rows/columns")
dev.off()
message("  saved -> 07_sample_distance_heatmap.pdf")


# ── 8. DESEQ2 RESULTS ─────────────────────────────────────────────────────────
# dds was already run with DESeq() in step 6
message("\n=== 8. DESeq2 results ===")

res_raw <- results(dds,
                   contrast = c("condition", "Cancer", "Normal"),
                   alpha    = PADJ_THRESHOLD)
summary(res_raw)

# ---- Shrinkage with apeglm --------------------------------------------------
# apeglm requires coef= not contrast=
# The coefficient name is built automatically from the design; check with:
#   resultsNames(dds)
coef_name <- grep("condition_Cancer", resultsNames(dds), value = TRUE)
if (length(coef_name) == 0) {
  message("  WARNING: Could not find 'condition_Cancer' coefficient.")
  message("  Available coefficients: ", paste(resultsNames(dds), collapse = ", "))
  message("  Using raw LFC (no apeglm shrinkage).")
  res_lfc <- res_raw
} else {
  res_lfc <- lfcShrink(dds, coef = coef_name[1], type = "apeglm")
}

# ---- Result table (base R,) ------------------------------------
res_df           <- as.data.frame(res_lfc)
res_df$gene      <- rownames(res_df)
res_df           <- res_df[order(res_df$padj, na.last = TRUE), ]
rownames(res_df) <- NULL   # clean up

res_df$DE_status <- ifelse(
  !is.na(res_df$padj) &
    res_df$padj < PADJ_THRESHOLD &
    res_df$log2FoldChange >= LFC_THRESHOLD,
  "UP in cancer",
  ifelse(
    !is.na(res_df$padj) &
      res_df$padj < PADJ_THRESHOLD &
      res_df$log2FoldChange <= -LFC_THRESHOLD,
    "DOWN in cancer",
    "not significant"
  )
)



cancer_samples <- colnames(dds)[dds$condition == "Cancer"]
normal_samples <- colnames(dds)[dds$condition == "Normal"]

# 1. Сырые counts (integer, до нормализации)
raw_mat <- counts(dds, normalized = FALSE)

mean_raw_cancer <- rowMeans(raw_mat[, cancer_samples, drop = FALSE])
mean_raw_normal <- rowMeans(raw_mat[, normal_samples, drop = FALSE])

# 2. Нормализованные counts (делённые на size factors)
norm_mat_all <- counts(dds, normalized = TRUE)

mean_norm_cancer <- rowMeans(norm_mat_all[, cancer_samples, drop = FALSE])
mean_norm_normal <- rowMeans(norm_mat_all[, normal_samples, drop = FALSE])

# 3. VST counts
vst_mat_all <- assay(vsd)

mean_vst_cancer <- rowMeans(vst_mat_all[, cancer_samples, drop = FALSE])
mean_vst_normal <- rowMeans(vst_mat_all[, normal_samples, drop = FALSE])


res_df$mean_raw_cancer  <- round(mean_raw_cancer[res_df$gene],  2)
res_df$mean_raw_normal  <- round(mean_raw_normal[res_df$gene],  2)

res_df$mean_norm_cancer <- round(mean_norm_cancer[res_df$gene], 2)
res_df$mean_norm_normal <- round(mean_norm_normal[res_df$gene], 2)

res_df$mean_vst_cancer  <- round(mean_vst_cancer[res_df$gene],  3)
res_df$mean_vst_normal  <- round(mean_vst_normal[res_df$gene],  3)


write.csv(res_df,
          file = file.path(OUTPUT_DIR, "DESeq2_results.csv"),
          row.names = FALSE)

n_up   <- sum(res_df$DE_status == "UP in cancer",   na.rm = TRUE)
n_down <- sum(res_df$DE_status == "DOWN in cancer",  na.rm = TRUE)
message("  UP in cancer: ", n_up, "   DOWN in cancer: ", n_down)


# ── 9. VISUALISATIONS ─────────────────────────────────────────────────────────
message("\n=== 9. Visualisations ===")

DE_COLORS <- c("UP in cancer"    = "#E53935",
               "DOWN in cancer"  = "#1E88E5",
               "not significant" = "grey72")

## 9.1  MA-plots (raw + shrunken, side by side) --------------------------------
pdf(file.path(OUTPUT_DIR, "08_MA_raw_vs_shrunken.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))
plotMA(res_raw, ylim = c(-8, 8), alpha = PADJ_THRESHOLD,
       main = "MA-plot  |  Raw LFC\n(noisy at low counts)")
plotMA(res_lfc, ylim = c(-8, 8), alpha = PADJ_THRESHOLD,
       main = "MA-plot  |  Shrunken LFC (apeglm)\nLow-count LFCs pulled towards 0")
dev.off()
message("  saved -> 08_MA_raw_vs_shrunken.pdf")

sig_genes <- res_df$gene[res_df$DE_status != "not significant"]
n_sig     <- length(sig_genes)
message("  Significant genes: ", n_sig)

if (n_sig >= 2) {
  top_genes_plot <- sig_genes
  plot_note      <- paste(n_sig, "significant DE genes")
} else {
  top_genes_plot <- head(res_df$gene, TOP_N)
  plot_note      <- paste("Top", TOP_N, "by padj (only", n_sig, "significant)")
  message("  Using top ", TOP_N, " by padj for all plots")
}


## 9.2  Volcano plot -----------------------------------------------------------
vol_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
vol_df$neg_log10_padj <- -log10(vol_df$padj + 1e-300)

top_vol <- vol_df[vol_df$gene %in% top_genes_plot, ]   

volcano_plot <- ggplot(vol_df, aes(x = log2FoldChange,
                                   y = neg_log10_padj,
                                   colour = DE_status)) +
  geom_point(alpha = 0.55, size = 1.3) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD),
             linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = top_vol, aes(label = gene),
                  size = 2.7, max.overlaps = 25) +
  scale_colour_manual(values = DE_COLORS) +
  labs(title    = "Volcano Plot  –  Cancer vs. Normal",
       subtitle = plot_note,    
       x        = "log2 Fold Change  (Cancer / Normal)",
       y        = "-log10(adjusted p-value)",
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
save_gg(volcano_plot, "09_volcano_plot", w = 12, h = 10)



## 9.3  Top-N individual gene expression (ggplot2 stripplot) -------------------
top_n_genes  <- top_genes_plot    

norm_mat     <- counts(dds, normalized = TRUE)[top_n_genes, , drop = FALSE]
norm_long_df <- data.frame(
  gene       = rep(rownames(norm_mat), times = ncol(norm_mat)),
  sample     = rep(colnames(norm_mat), each  = nrow(norm_mat)),
  norm_count = as.vector(norm_mat)
)
norm_long_df$Condition <- meta_data_filtered[norm_long_df$sample, "condition"]
norm_long_df$log2n     <- log2(norm_long_df$norm_count + 1)
norm_long_df$gene      <- factor(norm_long_df$gene, levels = rev(top_n_genes))

strip_plot <- ggplot(norm_long_df,
                     aes(x = log2n, y = gene, colour = Condition)) +
  geom_jitter(height = 0.2, size = 2.5, alpha = 0.8) +
  stat_summary(aes(group = Condition), fun = median,
               geom = "crossbar", width = 0.45, fatten = 2,
               colour = "black", linewidth = 0.35) +
  scale_colour_manual(values = c(Normal = "#43A047", Cancer = "#E53935")) +
  labs(title  = paste("Top", length(top_n_genes), "DE features – Normalised Expression"),
       subtitle = plot_note,    
       x      = "log2(normalised count + 1)",
       y      = NULL,
       colour = "Condition") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
save_gg(strip_plot, "10_top_genes_stripplot",
        h = max(5, length(top_n_genes) * 0.42))

## 9.4  Heatmap of 20 or less DE genes ------------------------------------------------
top_hm <- head(top_genes_plot, 20)

vst_counts <- assay(vsd)[top_hm, , drop = FALSE]

annotation_col <- data.frame(
  Condition = meta_data_filtered$condition,
  row.names = rownames(meta_data_filtered)
)
ann_colors <- list(Condition = c(Normal = "#43A047", Cancer = "#E53935"))

pdf(file.path(OUTPUT_DIR, "11_heatmap_top_DE.pdf"),
    width = 12, height = max(6, length(top_hm) * 0.4))
pheatmap(vst_counts,
         annotation_col    = annotation_col,
         annotation_colors = ann_colors,
         scale             = "row",
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 8,
         border_color      = NA,
         main              = paste("Top", length(top_hm),
                                   "DE features  |  VST + row Z-score"))
dev.off()



## 9.5  Top gene count plots ---------------------------------------------------


sig_only <- res_df[res_df$DE_status != "not significant", ]

if (nrow(sig_only) >= 2) {
 
  top_up   <- sig_only[which.max(sig_only$log2FoldChange), "gene"]
  top_down <- sig_only[which.min(sig_only$log2FoldChange), "gene"]
  rest     <- setdiff(sig_only$gene, c(top_up, top_down))
  rest     <- head(rest, 18)
  genes_95 <- c(top_up, top_down, rest)
  plot_note_95 <- paste0("top up (", top_up, ") + top down (", top_down, 
                         ") + ", length(rest), " significant")
  
} else {
 
  genes_95     <- head(res_df$gene, TOP_N)
  top_up       <- genes_95[1]
  top_down     <- genes_95[2]
  plot_note_95 <- paste0("top 20 by padj (only ", nrow(sig_only), " significant)")
  message("  9.5: < 2 significant genes, using top ", TOP_N, " by padj")
}

message("  9.5 genes: ", paste(genes_95, collapse = ", "))

# normalised comts (same as in stripplot)
norm_mat_95  <- counts(dds, normalized = TRUE)[genes_95, , drop = FALSE]
norm_long_95 <- data.frame(
  gene       = rep(rownames(norm_mat_95), times = ncol(norm_mat_95)),
  sample     = rep(colnames(norm_mat_95), each  = nrow(norm_mat_95)),
  norm_count = as.vector(norm_mat_95)
)
norm_long_95$Condition <- meta_data_filtered[norm_long_95$sample, "condition"]
norm_long_95$log2n     <- log2(norm_long_95$norm_count + 1)

# ORDE: top_up first , top_down second, rest are for padj
gene_order <- c(top_up, top_down, setdiff(genes_95, c(top_up, top_down)))
norm_long_95$gene <- factor(norm_long_95$gene, levels = rev(gene_order))

label_colors <- setNames(
  ifelse(gene_order == top_up,   "#E53935",
         ifelse(gene_order == top_down, "#1E88E5", "grey30")),
  gene_order
)

gene_plot_95 <- ggplot(norm_long_95,
                       aes(x = Condition, y = log2n, colour = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.6) +
  geom_jitter(position = position_jitter(width = 0.15), 
              size = 2, alpha = 0.6) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  scale_colour_manual(values = c(Normal = "#43A047", Cancer = "#E53935")) +
  labs(title    = "Top UP + DOWN + significant genes",
       subtitle = plot_note_95,
       x        = NULL,
       y        = "log2(normalised count + 1)",
       colour   = "Condition") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8, face = "bold"))

save_gg(gene_plot_95, "12_top_gene_counts",
        w = n_cols * 3.2,
        h = max(5, n_rows * 3.2))




# ── 10. SAVE WORKSPACE ────────────────────────────────────────────────────────
message("\n=== 10. Saving workspace ===")

save(dds, vsd, res_raw, res_lfc, res_df,
     volcano_plot, strip_plot, gene_plot_95,
     file = file.path(OUTPUT_DIR, "DESeq2_workspace.RData"))

message("\n✓  All done. Results are in: ", OUTPUT_DIR)
message(paste(" ", list.files(OUTPUT_DIR), collapse = "\n"))
