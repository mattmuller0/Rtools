###########################################################################
#
#                            expr_filtering_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-09-06
# Script Name: expr_filtering_functions

#======================== LIBRARIES ========================#
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(SummarizedExperiment)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')

#======================== CODE ========================#
# Function to get percent of genes detected
# Arguments:
#   dds: DESeq2 object
# Outputs:
#   data frame of percent of genes detected
percentGenesDetected <- function(dds, min_value = 0) {
  require(DESeq2)
  require(SummarizedExperiment)
  require(tidyverse)
  require(ggplot2)
  
  counts <- assay(dds)
  percent_genes_detected <- rowMeans(counts > min_value)
  return(percent_genes_detected)
}

# Function to run general preprocessing on a DESeq2 object
# Arguments:
#   dds: DESeq2 object
#   outpath: path to output directory
#   depth_filter: minimum library depth to keep
#   size_filter: minimum library size to keep
#   percent_filter: minimum percent of genes detected to keep
# Outputs:
#   plot of library depth
rna_preprocessing <- function(
  dds, outpath, 
  depth_filter = 10, size_filter = 1e6, percent_filter = 0.5, 
  normalize="none", group = NULL
  ) {
  require(DESeq2)
  require(SummarizedExperiment)
  require(tidyverse)
  require(ggplot2)
  require(edgeR)
  
  dir.create(outpath, showWarnings = F, recursive = T)
  # set normalization
  normal <- function(x) {normalize_counts(x, method = normalize)}

  # make a dds before filtering
  dds_before <- dds

  # Take a look at library depth
  message("Checking library depth")
  if (!is.na(depth_filter)) {
    before_plot <- plot_library_depth(dds_before, "Library Depth Prefiltering", bins = 100)
    keep <- rowMeans(as.data.frame(normal(dds))) >= depth_filter
    dds <- dds[keep,]
    after_plot <- plot_library_depth(dds, "Library Depth Postfiltering", bins = 100)
    depth_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  }
  if (is.na(depth_filter)) {
    before_plot <- plot_library_depth(dds, "Library Depth Prefiltering", bins = 100)
  }
  
  # Take a look at library size
  message("Checking library size")
  if (!is.na(size_filter)) {
    before_plot <- plot_library_size(dds_before, "Library Size Prefiltering", bins = 10)
    keep <- colSums(as.data.frame(normal(dds))) >= size_filter
    dds <- dds[,keep]
    after_plot <- plot_library_size(dds, "Library Size Postfiltering", bins = 10)
    size_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  }
  if (is.na(size_filter)) {
    before_plot <- plot_library_size(dds, "Library Size Prefiltering", bins = 10)
  }

  # Take a look at percent of genes detected
  message("Checking percent of genes detected")
  if (!is.na(percent_filter)) {
    before_plot <- plot_percent_genes_detected(dds_before, "Percent of Genes Detected Prefiltering")
    keep <- percentGenesDetected(dds) >= percent_filter
    dds <- dds[keep,]
    after_plot <- plot_percent_genes_detected(dds, "Percent of Genes Detected Postfiltering")
    percent_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  }
  if (is.na(percent_filter)) {
    before_plot <- plot_percent_genes_detected(dds, "Percent of Genes Detected Prefiltering")
  }

  # save all the plots together
  all_plots <- plot_grid(depth_plots, size_plots, percent_plots, nrow = 3)
  ggsave(file.path(outpath, "preprocessing_plots.pdf"), all_plots)
  
  # Make a tree
  message("Checking hierarchical clustering")
  counts <- scale(t(normal(dds)))
  sampleTree = hclust(dist(counts))
  ggtree(sampleTree) + 
    geom_tiplab(size = 2, hjust = -0.1) + 
    theme_tree2() +
    labs(title = "Sample clustering to detect outliers", subtitle = "", x = "", y = "")
  ggsave(file.path(outpath, "sample_outliers.pdf"))

  # PCA Plot
  message("Checking PCA")
  counts <- scale(t(assay(vst(dds))))
  pca <- prcomp(counts)
  if (is.null(group)) {pca_plot <- ggbiplot(pca, ellipse = TRUE, var.axes = FALSE)}
  if (!is.null(group)) {pca_plot <- ggbiplot(pca, groups = colData(dds)[, group], ellipse = TRUE, var.axes = FALSE)}

  pca_plot <- pca_plot +
    theme(legend.direction = 'horizontal', legend.position = 'top') +
    geom_text_repel(aes(label = rownames(counts)), size = 4, vjust = -1, max.overlaps = log10(length(colnames(dds)))) +
    theme_classic2() +
    labs(title = "PCA Plot")
  ggsave(file.path(outpath, "pca_plot.pdf"), pca_plot, dpi = 300)

  # save the dds object
  message("Saving dds object")
  saveRDS(dds, file = file.path(outpath,"dds.rds"))
  
  # we want at least the raw and vst normalized counts
  save_se(dds, outpath, normalize = normalize)
  if (normalize != "none") {
    save_se(dds, outpath, normalize = "none")
  }
  if (normalize != "vst") {
    save_se(dds, outpath, normalize = "vst")
  }

  return(dds)
}

# Function to filter a dds by expression and make plots of the filtering
# Arguments:
#   dds: DESeq2 object
#   outpath: path to output directory
#   group: column of interest
#   min.count: minimum mean expression to keep
#   min.prop: minimum proportion of samples to keep
#   ...: additional arguments to pass to FilterByExpr
filter_edgeR <- function(
  dds, outpath, 
  group = NULL, 
  min.count = 10, min.prop = 0.5,
  ...
  ) {
  require(DESeq2)
  require(SummarizedExperiment)
  require(tidyverse)
  require(ggplot2)
  require(edgeR)

  dir.create(outpath, showWarnings = F, recursive = T)

  # make a dds before filtering
  dds_before <- dds

  # filter by expression
  message("Filtering with FilterByExpr from edgeR")
  keep <- filterByExpr(dds, group = dds[[group]], min.count = min.count, ...)
  dds <- dds[keep,]

  # Plot the filtering
  message("Plotting filtering results")
  # Plot the library depth prefiltering
  before_plot <- plot_library_depth(dds_before, "Library Depth Prefiltering", bins = 100)
  # Plot the library depth postfiltering
  after_plot <- plot_library_depth(dds, "Library Depth Postfiltering", bins = 100)

  # Plot the library size prefiltering
  before_plot2 <- plot_library_size(dds_before, "Library Size Prefiltering", bins = 10)
  # Plot the library size postfiltering
  after_plot2 <- plot_library_size(dds, "Library Size Postfiltering", bins = 10)

  # Plot the percent of genes detected prefiltering
  before_plot3 <- plot_percent_genes_detected(dds_before, "Percent of Genes Detected Prefiltering")
  # Plot the percent of genes detected postfiltering
  after_plot3 <- plot_percent_genes_detected(dds, "Percent of Genes Detected Postfiltering")

  # Plot all the plots together
  depth_plots <- plot_grid(before_plot, after_plot, ncol = 2)
  size_plots <- plot_grid(before_plot2, after_plot2, ncol = 2)
  percent_plots <- plot_grid(before_plot3, after_plot3, ncol = 2)
  all_plots <- plot_grid(depth_plots, size_plots, percent_plots, nrow = 3)
  # Save the plots
  ggsave(file.path(outpath, "filtering_plots.pdf"), all_plots)

  # Make a tree
  message("Checking hierarchical clustering")
  counts <- scale(t(assay(vst(dds))))
  sampleTree = hclust(dist(counts))
  ggtree(sampleTree) + 
    geom_tiplab(size = 2, hjust = -0.1) + 
    theme_tree2() +
    labs(title = "Sample clustering to detect outliers", subtitle = "", x = "", y = "")
  ggsave(file.path(outpath, "sample_outliers.pdf"), dpi = 300)

  # PCA Plot
  message("Checking PCA")
  # get the variance stabilized counts
  counts <- scale(t(assay(vst(dds))))
  pca <- prcomp(counts)
  if (is.null(group)) {pca_plot <- ggbiplot(pca, ellipse = TRUE, var.axes = FALSE)}
  if (!is.null(group)) {pca_plot <- ggbiplot(pca, groups = colData(dds)[, group], ellipse = TRUE, var.axes = FALSE)}
  # plot the pca
  pca_plot <- pca_plot +
    theme(legend.direction = 'horizontal', legend.position = 'top') +
    geom_text_repel(aes(label = rownames(counts)), size = 4, vjust = -1, max.overlaps = log10(length(colnames(dds)))) +
    theme_classic2() +
    labs(title = "PCA Plot")
  ggsave(file.path(outpath, "pca_plot.pdf"), pca_plot, dpi = 300)

  # save the dds object
  message("Saving dds object")
  saveRDS(dds, file = file.path(outpath,"dds.rds"))

  # we want at least the raw and vst normalized counts
  save_se(dds, outpath, normalize = "none")
  save_se(dds, outpath, normalize = "vst")

  return(dds)
}

# Function to detect wbc contamination in platelet rna-seq
# Arguments:
#   dds: DESeq2 object
#   outpath: path to output directory
#   group: column of interest
# Returns:
#   dds: DESeq2 object
#   plot of wbc contamination
detect_wbc <- function(dds, outpath, group = NULL, normalize = "mor") {
  require(DESeq2)
  require(SummarizedExperiment)
  require(tidyverse)
  require(ggplot2)
  require(edgeR)

  dir.create(outpath, showWarnings = F, recursive = T)

  # make a dds before filtering
  dds_before <- dds

  # get the counts
  counts <- normalize_counts(dds, normalize)

  # get the ratio
  ratio <- counts['PTPRC',] / counts['ITGA2B',]

  # get the group
  tryCatch(
    grouping <- colData(dds)[,group],
    error = function(e) {
      message("No group specified or group does not exist")
      grouping <- ratio > 1
    }
  )
  IDs <- rownames(colData(dds))

  plotting_df <- data.frame(IDs, ratio, grouping)

  # plot the ratio
  plot <- ggplot(data.frame(ratio), aes(x = IDs, y = ratio, col = grouping)) +
      geom_point() +
      theme_matt(18) +
      labs(title = "WBC Contamination", x = "PTPRC/ITGA2B", y = "Count") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      geom_text_repel(label = IDs, size = 5, show.legend = FALSE)

  ggsave(file.path(outpath, "wbc_contamination.pdf"), plot)
  write.csv(plotting_df, file.path(outpath, "wbc_contamination.csv"))
  saveRDS(dds, file.path(outpath,"dds.rds"))

  return(dds)
}