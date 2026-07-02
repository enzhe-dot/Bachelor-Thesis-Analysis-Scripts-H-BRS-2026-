# ============================================================
# ORA (Over-Representation Analysis) for targets miRNA
# 4 cancer types: Breast, Lung, Larynx, Hepatocellular
# GO (BP, MF, CC) + KEGG
# Visualisation: dotplot, emapplot, cnetplot
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# ---- 1. Paths ----
output_dir <- "C:/Users/ПК/OneDrive/Документы/ORA_targets_v2"

#adjust ofc 
cancer_files <- c(
  Breast         = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Breast.csv",
  Lung           = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Lung.csv",
  Larynx         = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Larynx.csv",
  Hepatocellular = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Hepatocellular.csv"
)


plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)





# ---- 3. Functions of ORA for target  ----
run_ORA <- function(cancer_name, file_path, universe) {
  
  cat(sprintf("\n========== %s ==========\n", cancer_name))
  
  
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Unique Entrez ID
  genes <- as.character(unique(df$Entrez_ID))
  genes <- genes[!is.na(genes) & genes != ""]
  cat(sprintf("Gene-targets: %d\n", length(genes)))
  
  results <- list()
  
  # ---- GO: Biological Process ----
  cat("  GO BP...\n")
  ego_BP <- tryCatch(
    enrichGO(gene          = genes,
             OrgDb         = org.Hs.eg.db,
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             minGSSize     = 10,
             maxGSSize     = 500,
             readable      = TRUE),   
    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
  )
  results$GO_BP <- ego_BP
  
  # ---- GO: Molecular Function ----
  cat("  GO MF...\n")
  ego_MF <- tryCatch(
    enrichGO(gene          = genes,
             OrgDb         = org.Hs.eg.db,
             ont           = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             minGSSize     = 10,
             maxGSSize     = 500,
             readable      = TRUE),
    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
  )
  results$GO_MF <- ego_MF
  
  # ---- GO: Cellular Component ----
  cat("  GO CC...\n")
  ego_CC <- tryCatch(
    enrichGO(gene          = genes,
             #universe      = universe,
             OrgDb         = org.Hs.eg.db,
             ont           = "CC",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             minGSSize     = 10,
             maxGSSize     = 500,
             readable      = TRUE),
    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
  )
  results$GO_CC <- ego_CC
  
  # ---- KEGG ----
  cat("  KEGG...\n")
  ekegg <- tryCatch(
    enrichKEGG(gene          = genes,
               organism      = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               minGSSize     = 10,
               maxGSSize     = 500),
    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
  )
  
  if (!is.null(ekegg)) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  results$KEGG <- ekegg
  
  return(results)
}


# ---- 4. Save tables  ----
save_tables <- function(ora_results, cancer_name, tables_dir) {
  
  for (analysis_name in names(ora_results)) {
    res <- ora_results[[analysis_name]]
    if (is.null(res)) next
    
    df <- as.data.frame(res)
    if (nrow(df) == 0) {
      cat(sprintf("  [%s %s] no \n", cancer_name, analysis_name))
      next
    }
    
    # Sorting p.adjust
    df <- df[order(df$p.adjust), ]
    
    out_path <- file.path(tables_dir, 
                          sprintf("%s_%s_ORA.csv", cancer_name, analysis_name))
    write.csv(df, out_path, row.names = FALSE)
    cat(sprintf("  [%s %s] %d significant pathways → %s\n", 
                cancer_name, analysis_name, nrow(df), basename(out_path)))
  }
}


# ---- 5.Visualisation function----
save_plots <- function(ora_results, cancer_name, plots_dir) {
  
  for (analysis_name in names(ora_results)) {
    res <- ora_results[[analysis_name]]
    if (is.null(res)) next
    
    df <- as.data.frame(res)
    if (nrow(df) == 0) next
    
    n_show <- min(20, nrow(df)) 
    
    # -- Dotplot --
    tryCatch({
      p <- dotplot(res, showCategory = n_show, 
                   title = sprintf("%s | %s | Dotplot (top %d)", 
                                   cancer_name, analysis_name, n_show)) +
        theme(axis.text.y = element_text(size = 8))
      
      ggsave(file.path(plots_dir, 
                       sprintf("%s_%s_dotplot.png", cancer_name, analysis_name)),
             plot = p, width = 10, height = 8, dpi = 150)
    }, error = function(e) cat(sprintf("  dotplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
    
    # -- Emapplot 
    if (nrow(df) >= 2) {
      tryCatch({
        res2 <- pairwise_termsim(res)
        p2 <- emapplot(res2, showCategory = n_show,
                       title = sprintf("%s | %s | Network", cancer_name, analysis_name))
        ggsave(file.path(plots_dir,
                         sprintf("%s_%s_emapplot.png", cancer_name, analysis_name)),
               plot = p2, width = 10, height = 9, dpi = 150)
      }, error = function(e) cat(sprintf("  emapplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
    }
    
    # -- Cnetplot — top 5 pathways --
    tryCatch({
      p3 <- cnetplot(res, showCategory = 5,
                     title = sprintf("%s | %s | Gene-Concept Network", 
                                     cancer_name, analysis_name))
      ggsave(file.path(plots_dir,
                       sprintf("%s_%s_cnetplot.png", cancer_name, analysis_name)),
             plot = p3, width = 12, height = 10, dpi = 150)
    }, error = function(e) cat(sprintf("  cnetplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
    
    # -- Barplot (alternative) --
    tryCatch({
      p4 <- barplot(res, showCategory = n_show,
                    title = sprintf("%s | %s | Barplot", cancer_name, analysis_name)) +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(plots_dir,
                       sprintf("%s_%s_barplot.png", cancer_name, analysis_name)),
             plot = p4, width = 10, height = 8, dpi = 150)
    }, error = function(e) cat(sprintf("  barplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
  }
}


# ---- 6.Launching ----
all_results <- list()

for (cn in names(cancer_files)) {
  all_results[[cn]] <- run_ORA(cn, cancer_files[cn], universe_genes)
  save_tables(all_results[[cn]], cn, tables_dir)
  save_plots(all_results[[cn]], cn, plots_dir)
}


# ---- 7.comparative dotplot — all cancers together (только GO BP) ----


cat("\n Comparative plot all cancers (GO BP)...\n")

tryCatch({
  bp_list <- lapply(names(all_results), function(cn) {
    all_results[[cn]][["GO_BP"]]
  })
  names(bp_list) <- names(all_results)
  
  
  bp_list <- bp_list[!sapply(bp_list, is.null)]
  
  if (length(bp_list) >= 2) {
    merged <- merge_result(bp_list)
    
    p_compare <- dotplot(merged, showCategory = 20,
                         title = "GO BP — comparison on cancers") +
      theme(axis.text.y = element_text(size = 8))
    
    ggsave(file.path(plots_dir, "ALL_cancers_GOBP_comparison_dotplot.png"),
           plot = p_compare, width = 14, height = 10, dpi = 150)
    cat("  Saved comparative dotplot\n")
  }
}, error = function(e) cat("  comparative plot:", e$message, "\n"))


# ---- 8. Sum up ----
cat("\n===== ENDing =====\n")
for (cn in names(all_results)) {
  cat(sprintf("\n[%s]\n", cn))
  for (an in names(all_results[[cn]])) {
    res <- all_results[[cn]][[an]]
    if (!is.null(res)) {
      n <- nrow(as.data.frame(res))
      cat(sprintf("  %s: %d significant pathways\n", an, n))
    } else {
      cat(sprintf("  %s: no results\n", an))
    }
  }
}

cat(sprintf("\nTables: %s\n", tables_dir))
cat(sprintf("Graphs: %s\n", plots_dir))