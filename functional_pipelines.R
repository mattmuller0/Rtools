###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-18
# 
# Script Name: Functional Pipelines
# 
# Notes:
# This script contains functions for running pipelines within an R interface. These functions mainly relate to rna-sequencing, but this is not exclusive.


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "devtools",
  "BiocManager",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "AnnotationDbi",
  "apeglm",
  "ggpubr",
  "EnhancedVolcano",
  "cowplot",
  "ggtree"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T) # nolint
}

# LOAD FUNCTIONS
# space reserved for sourcing in functions
options(stringsAsFactors = FALSE)

source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/stats_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/deseq_helper_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/geneset_creation_functions.R')

###########################################################################
#
#                                 CODE
#
###########################################################################
##### Geneset Creation Pipeline
# Input:
#   - outdir: output directory
#   - dds: DESeq2 object
#   - target: target variable
#   - depthFilter: minimum number of samples a gene must be expressed in
#   - sizeFilter: minimum number of genes a geneset must contain
#   - percentFilter: minimum percent of genes in a geneset that must be expressed in a sample
#   - test_names: names of tests to run for geneset creation
#   - test_params: parameters for each test
create_geneset <- function(
  # metadata
  outdir,
  # data 
  dds, target, 
  # rna filtering parameters
  depthFilter=10, sizeFilter=1e6, percentFilter=0.25,
  # geneset creation parameters 
  deseq_filter = 0.05, fc_filter = 0.5,
  cor_filter = NA, cor_method = "spearman",
  # predictive modeling parameters
  normalization = "cpm",
  # sys params
  conda = "tf", kwargs = NULL
  ) {
  # rna preprocessing
  message("========= Preprocessing Data =========")
  message("Filtering Data...")
  dds <- rna_preprocessing(dds, file.path(outdir, "preprocessing"), depthFilter, sizeFilter, percentFilter)

  # DESeq analysis
  message("========= Running DESeq Analysis =========")
  message(paste0("Pulling DESeq results for: ", target, "..."))
  threshold <- "padj"
  deseq_res_name <- c(target, "1", "0")
  deseq_res <- results(dds, contrast = deseq_res_name) %>% 
    as.data.frame()
  deseq_sign_genes <- deseq_res %>%
    filter(padj < deseq_filter, abs(log2FoldChange) > fc_filter) %>% 
    rownames()
  # If not enough genes identified by FDR correction, switch to pvalue threshold
  if (length(deseq_sign_genes) < 75) {
    message("Not enough genes identified by FDR correction. Switching to pvalue threshold.")
    deseq_sign_genes <- deseq_res %>% 
      filter(pvalue < deseq_filter, abs(log2FoldChange) > fc_filter) %>% 
      rownames()
      threshold <- "pvalue"
      }
  message(paste0("DESeq identified genes: n = ", length(deseq_sign_genes)))
  # Correlation analysis
  message("========= Running Correlation Analysis =========")
  message(paste0("Pulling Correlation results for: ", target, "..."))
  cor_res <- calculate_correlations(dds, target, method=cor_method)
  cor_sign_genes <- cor_res %>% 
    filter(abs(correlation) > cor_filter) %>% 
    rownames()
  message(paste0("Correlation identified genes: n = ", length(cor_sign_genes)))

  # Create geneset
  message("========= Creating Geneset =========")
  message(paste0("Creating geneset for: ", target, "..."))
  geneset_genes <- intersect(deseq_sign_genes, cor_sign_genes)
  message(paste0("Geneset contains: ", length(geneset_genes), " genes"))
    res <- deseq_res %>% dplyr::select(baseMean, log2FoldChange, pvalue, padj) %>%
      mutate(
        deseq_pval = pvalue, deseq_padj = padj,
        correlation = cor_res$correlation, cor_pval = cor_res$pvalue, cor_padj = cor_res$padj,
        pvalue = NULL, padj = NULL
        )

  # Plot metrics
  message("========= Plotting Results =========")
  message("Plotting geneset volcano...")
  max_log2FC <- max(res$log2FoldChange)
  min_log2FC <- min(res$log2FoldChange)
  geneset_correlations <- EnhancedVolcano(
    res, lab=rownames(res),
    x = "log2FoldChange", y = threshold,
    pCutoff = deseq_filter, FCcutoff = fc_filter
    ) + 
    labs(title = paste0("Geneset Volcano: ", target), x = "log2FoldChange", y = threshold)
  ggsave(file.path(outdir, paste0("geneset_volcano_", target, ".png")), geneset_correlations)

  message("Plotting geneset correlations...")
  geneset_correlations <- EnhancedVolcano(
    res, lab=rownames(res),
    x = "log2FoldChange", y = "correlation",
    pCutoff = cor_filter, FCcutoff = fc_filter
    ) +
    lims(y = c(-1, 1)) +
    labs(title = paste0("Geneset Correlations: ", target), x = "log2FoldChange", y = "correlation")
    ggsave(file.path(outdir, paste0("geneset_correlations_", target, ".png")), geneset_correlations)

  # export geneset
  message("Exporting geneset...")
  junk <- export_counts_and_labels_from_SE(dds, geneset_genes, target, outdir, normalize = normalization)

  # feature elimination and predictive modeling
  message("========= Running Predictive Modeling =========")
  sys.call(paste0("conda activate ", conda))
  sys.call("Python3", c(
    "recursive_feature_elimination_and_model_generation.py",
    "--input", file.path(outdir, paste0(target, "_counttable")),
    "--labels", file.path(outdir, paste0(target, "_labels")),
    "--output", file.path(outdir, "modeling"),
    "--scoring", "roc_auc",
    "--random_state", 420,
    "--folds", 5,
    "--verbose", TRUE
    ))
  sys.call("conda deactivate")

  # Load results
  message("========= Loading Predictive Modeling Results =========")
  message("Loading results...")
  # Load in the final geneset
  geneset_final <- read.csv(paste0(outdir, "/modeling/optimal_features.csv"))$X0
  # Load in the test and train cohorts
  test <- read.csv(paste0(outdir, "/modeling/data/X_test.csv"))
  train <- read.csv(paste0(outdir, "/modeling/data/X_train.csv"))

  # Stop early if geneset is too small
  if (length(geneset_final) < 10) {
    message("Geneset is too small. Stopping early.")
    return()
  }

  IRkernel::installspec()
  message("Plotting geneset heatmap...")
  geneset_heatmap <- plot_geneset_heatmap(dds, geneset_final, target)
  ggsave(file.path(outdir, paste0("geneset_heatmap_", target, ".png")), geneset_heatmap)

  message("Plotting geneset PCA...")
  geneset_pca <- plot_normalized_pca(dds[geneset_final,], target, normalization)
  ggsave(file.path(outdir, paste0("geneset_pca_", target, ".png")), geneset_pca)

  message("Plotting geneset boxplots...")
  geneset_boxplots <- plot_gene_boxplot(dds, geneset_final, target, file.path(outdir, "geneset_boxplots"))
  ggsave(file.path(outdir, paste0("geneset_boxplots_", target, ".png")), geneset_boxplots)

  message("Plotting geneset singscores...")
  up_genes <- hyper_res %>% 
    filter(log2FoldChange > 0) %>%
    select(geneset_final) %>%
    row.names()
  down_genes <- hyper_res %>%
    filter(log2FoldChange < 0) %>%
    select(geneset_final) %>%
    row.names()
  geneset_singscores <- plot_singscores(dds[,rownames(test)], target, file.path(outdir, "geneset_singscores.png"), up_genes, down_genes)

  message("Saving parameters...")
  params <- data.frame(
    target = target,
    geneset = geneset,
    depthFilter = depthFilter,
    sizeFilter = sizeFilter,
    percentFilter = percentFilter,
    deseq_filter = deseq_filter,
    fc_filter = fc_filter,
    cor_filter = cor_filter,
    cor_method = cor_method,
    normalization = normalization,
    conda_env = conda
    )
  write.csv(params, file.path(outdir, paste0(target, "_modeling_params.csv")), row.names = F)
  }
