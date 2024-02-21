###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-03-10
# 
# Script Name: Geneset creation Functions
# 
# Notes:
# This is a set of functions used to create gene sets that can be used in
# a variety of gene set enrichment analyses. The functions are designed to
# be used with SummarizedExperiment objects, but can be used with other
# data types as well.

# TODO: Add in support for clusterprofiler gene set enrichment analysis

###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
library(tidyverse)
library(ggplot2)
library(SummarizedExperiment)
library(DESeq2)
library(singscore)
library(ggbiplot)
library(edgeR)
library(EnhancedVolcano)
library(ggpubr)
library(GGally)
library(ggrepel)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/deseq_helper_functions.R')

###########################################################################
#
#                                 CODE
#
###########################################################################
# Function to create a gene set from a DESeq2 object
# Inputs:
#   dds: DESeq2 object
#   binary_condition: name of the binary condition to use
#   outpath: path to output directory
#   deseq_filter: p-value filter for DESeq2 results
#   wilcox_filter: p-value filter for Wilcoxan results
#   cor_filter: p-value filter for correlation results
#   fc_filter: fold change filter for DESeq2 results
#   deseq_res_name: name of the DESeq2 results to use
#   correlation_method: method to use for correlation
#   correlation_var: variable name to use for correlations
#   corrected: whether or not to use the corrected p-value
# Outputs:
#   The genes and statistics for the geneset
geneset_creation <- function(
    # Inputs
    dds, binary_condition, outpath,
    # Filters
    deseq_filter = NULL, wilcox_filter = NULL, cor_filter = NULL, fc_filter = 0.5, 
    # Methodology
    deseq_res_name = c(binary_condition, "1", "0"), 
    correlation_method = "pearson", 
    correlation_var = NULL,
    corrected = F,

    # settings
    image_type = "pdf"
    ) {
  require(tidyverse)
  require(EnhancedVolcano)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  require(singscore)

  # make sure directory exists
  dir.create(outpath, showWarnings = F, recursive = T)

  # Make sure the correlation variable is set
  if (is.null(correlation_var)) {correlation_var <- binary_condition}
  
  # DESeq analysis
  message("Running DESeq Analysis")
  print(paste0("Pulling DESeq results for: ", binary_condition, "..."))
  deseq_res <- results(dds, contrast = deseq_res_name) %>% 
    as.data.frame()
  if (!corrected) {deseq_sign_genes <- deseq_res %>% filter(pvalue < deseq_filter, abs(log2FoldChange) > fc_filter) %>% rownames()}
  if (corrected) {deseq_sign_genes <- deseq_res %>% filter(padj < deseq_filter, abs(log2FoldChange) > fc_filter) %>% rownames()}
  message("Finished DESeq Analysis")
  print(paste0("DESeq identified genes: n = ", length(deseq_sign_genes)))

  # Wilcoxan Analysis
  message("Running Wilcox Analysis")
  wilcoxan_res <- gene_wilcox_test(dds, binary_condition)
  if (!corrected) {wilcox_sign_genes <- wilcoxan_res %>% filter(pvalue < wilcox_filter) %>% rownames}
  if (corrected) {wilcox_sign_genes <- wilcoxan_res %>% filter(padj < wilcox_filter) %>% rownames}  
  message("Finished Wilcox Analysis")
  print(paste0("Wilcox identified genes: n = ", length(wilcox_sign_genes)))


  # Correlation Analysis
  message("Running Correlation Analysis")
  cor_res <- calculate_correlations(dds, correlation_var, method=correlation_method)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = 'fdr')
  # cor_sign_genes <- cor_res %>% filter(abs(correlation) > 0.2) %>% rownames
  cor_sign_genes <- cor_res %>% 
    filter(pvalue < cor_filter) %>% 
    rownames()
  message("Finished Correlation Analysis")
  print(paste0("Correlation identified genes: n = ", length(cor_sign_genes)))
  
  message("Making Geneset Output")
  geneset_df <- deseq_res %>% dplyr::select(baseMean, log2FoldChange, pvalue, padj) %>%
    mutate(deseq_pval = pvalue, deseq_padj = padj,
          wilcox_pval = wilcoxan_res$pvalue, wilcox_padj = wilcoxan_res$padj,
          correlation = cor_res$correlation, cor_pval = cor_res$pvalue, cor_padj = cor_res$padj,
          pvalue = NULL, padj = NULL)

  # Make some metrics on the results
  cor_plot <- ggplot(geneset_df, aes(x=correlation, y=-log10(cor_pval))) + 
        geom_point(aes(col = cut(cor_pval, c(-Inf, 0.05, Inf)))) +
        scale_color_manual(name = "pvalue",
                            values = c("(-Inf,0.05]" = "red",
                                      "(0.05, Inf]" = "black"),
                            labels = c("pvalue <= 0.05", "NS")) +
        geom_text_repel(aes(label = ifelse(cor_pval <= 0.01, rownames(geneset_df), '')), size = 3, vjust = -1) +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        labs(title = "Expression Correlation Plot", y = '-log(pvalue)') + 
        theme_classic2()
  ggsave(paste0(outpath, "correlation_pval_plot.png"), cor_plot)

  # Correlation plot to logFC
  cor_plot <- ggplot(geneset_df, aes(x=log2FoldChange, y=correlation)) + 
        geom_point(aes(col = cut(cor_pval, c(-Inf, 0.05, Inf))), size = 1) +
        scale_color_manual(name = "pvalue",
                            values = c("(-Inf,0.05]" = "red",
                                      "(0.05, Inf]" = "black"),
                            labels = c("pvalue <= 0.05", "NS")) +
        geom_text_repel(aes(label = ifelse(cor_pval <= 0.01, rownames(geneset_df), '')), size = 3, vjust = -1) +
        labs(title = "Expression Correlation Plot", x = "Log[2] Fold Change") + 
        theme_classic2()
  ggsave(paste0(outpath, "correlation_logfc_plot.png"), cor_plot)

  # Plot Wilcoxan
  wilcox_plot <- EnhancedVolcano(geneset_df, lab = rownames(geneset_df),
                            x = 'log2FoldChange', y = 'wilcox_pval',
                            pCutoff = wilcox_filter, FCcutoff = fc_filter)
  ggsave(paste0(outpath, "wilcox_plot.", image_type), wilcox_plot)
  
  # Venn Diagram
  png(paste0(outpath, 'venn_diagram.png'))
  venn_list <- list(deseq = deseq_sign_genes, wilcox = wilcox_sign_genes, correlation = cor_sign_genes)
  # remove with 0 genes or NULL
  sets_to_use <- sapply(venn_list, function(x) length(x) > 0)
  venn_list <- venn_list[sets_to_use]
  vennD <- venn(venn_list)
  dev.off()

  # enrichment
  x <- deseq_res %>% dplyr::select(log2FoldChange) %>% arrange(-log2FoldChange)
  geneList <- x$log2FoldChange
  names(geneList) <- row.names(x)
  gse <- rna_enrichment(geneList, outpath,  keyType = "SYMBOL")

  # save deseq results
  write.csv(deseq_res, file = paste0(outpath, "deseq_results.csv"))

  # Get the prelim geneset using reduce
  prelim_geneset <- Reduce(intersect, venn_list)

  # Filter to genes and return
  geneset_df <- geneset_df %>% filter(row.names(geneset_df) %in% prelim_geneset)
  message(paste0("Geneset N genes: ", length(prelim_geneset) ) )
  message("All Done!")
  return(geneset_df)
}



###########################################################################
#
#                       Validation Functions
#
###########################################################################
# Function to plot the singscore across a condition
plot_singscore <- function(dds, condition, up_genes, down_genes, jitter=FALSE) {
  scoring <- simpleScore(
    dds %>% rankGenes, 
    upSet = up_genes, downSet = down_genes
    ) %>% 
      mutate(condition = colData(dds)[,condition])
  singscore_plot <- ggplot(scoring, aes(x = condition, y = TotalScore, fill = condition)) + 
    geom_boxplot() + 
    labs(title = "Singscore Boxplot", x = "Condition") + 
    theme_classic2()
  if (jitter) {singscore_plot <- singscore_plot + geom_jitter()}
  return(singscore_plot)
}

# Function to plot the gene counts across a condition in a boxplot of a given geneset
plot_gene_boxplot <- function(dds, geneset, condition, log_scaled = T) {
  require(reshape2)
  require(ggplot2)
  require(tidyverse)
  require(ggpubr)  

  cond <- colData(dds)[,condition]
  pace_counts <- dds[geneset,] %>% counts(normalize=T) %>% t() %>% as.data.frame()
  if (log_scaled) {pace_counts <- log2(pace_counts)}
  
  pace_counts$condition <- cond
  data <- pace_counts %>% melt() %>% arrange(condition, variable, value)
  
  plot <- ggplot(data, aes(x=variable, y=value, fill=condition)) + 
    geom_boxplot() + 
    labs(title = "Geneset Boxplot", y = "Counts", x = "Genes") + 
    theme_classic2() +
    theme(axis.text.x = element_text(angle = -45, hjust = 0))
  return(plot)
}

export_counts_and_labels_from_SE <- function(dds, geneset, label, outpath, normalize="mor") {
  require(DESeq2)
  require(edgeR)
  require(tidyverse)
  
  if (normalize == "mor") {counttable <- counts(dds, normalize = T)}
  if (normalize %in% c("vsd", "vst")) {counttable <- assay(varianceStabilizingTransformation(dds))}
  if (normalize == "log2") {counttable <- log2(assay(dds))}
  if (normalize == "rld") {counttable <- rlog(dds)}
  if (normalize == "cpm") {counttable <- cpm(dds)}
  
  
  counttable <- counttable[geneset,]
  labels <- colData(dds)[, label]
  names(labels) <- colnames(counttable)
  
  message(paste0("Writing output counttable and labels to ", outpath))
  write.csv(counttable, file = paste0(outpath, "/", label, "_counttable.csv"), row.names = T)
  write.csv(labels, file = paste0(outpath, "/", label, "_labels.csv"), row.names = T, quote = F)
  return(list(counttable, labels))
}

plot_normalized_pca <- function(dds, group, normalize="mor", labels = FALSE) {
  # Make a normalized PCA plot
  if (normalize == "mor") {counttable <- counts(dds, normalize = T)}
  if (normalize %in% c("vsd", "vst")) {counttable <- assay(varianceStabilizingTransformation(dds))}
  if (normalize == "log2") {counttable <- log2(assay(dds)+1)}
  if (normalize == "rld") {counttable <- rlog(dds)}
  if (normalize == "vst") {counttable <- rlog(dds)}
  if (normalize == "cpm") {counttable <- cpm(dds)}
  
  # get the group
  condition <- colData(dds)[,group]

  pca <- prcomp(scale(t(counttable)))
  pca_plot <- ggbiplot(pca, group = condition, var.axes = F) + 
    theme_classic2() + 
    labs(title = paste0("PCA Plot"))
  if (labels) {pca_plot <- pca_plot + geom_text_repel(aes(label = rownames(counttable)), size = 2)}
  return(pca_plot)
}

# Prep data to match the geneset for the validation
prep_validation_data <- function(dds, geneset, label, outpath, normalize="mor") {
  # Normalize the data
  counttable <- normalize_counts(dds, normalize)

  # Subset and prep the data
  counttable <- counttable %>% add_missing_rows(geneset)
  counttable <- counttable[geneset,]

  # Add the labels
  labels <- colData(dds)[, label]
  names(labels) <- colnames(counttable)

  # export the data for use in the geneset prediction pipeline
  write.csv(counttable, file = file.path(outpath, 'validation_', label, '_data.csv'))
  write.csv(labels, file = file.path(outpath, 'validation_', label, '_labels.csv'))

  # save the dds object
  saveRDS(dds, file = file.path(outpath, 'validation_', label, '_dds.rds'))
  }

plot_geneset_heatmap <- function(dds, geneset, label) {
  require(DESeq2)
  require(edgeR)
  require(tidyverse)
  require(gplots)
  require(ComplexHeatmap)
  
  # Subset and prep the data
  counttable <- counts(dds, normalize = T)
  counttable <- counttable %>% add_missing_rows(geneset)
  counttable <- counttable[geneset,]
  counttable <- log2(counttable+1)
  
  # Add the labels
  labels <- colData(dds)[, label]
  names(labels) <- colnames(counttable)
  
  # Make the heatmap
  hm <- Heatmap(counttable)
  return(hm)
}