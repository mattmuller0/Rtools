###########################################################################
#
#                            rna_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: rna_functions

#======================== LIBRARIES ========================#
library(lintr) #nolint
library(httpgd) #nolint
library(languageserver) #nolint
library(devtools)
library(tidyverse)
library(glue)
library(AnnotationDbi)
library(enrichplot)
library(forcats)
library(ggplot2)
library(msigdbr)
library(ggpubr)
library(ggbiplot)
library(ggpubr)
library(DESeq2)
library(edgeR)
library(limma)


# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/plotting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/enrichment_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/converting_functions.R')

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
  
  counts <- dds %>% assay()
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
    # ggsave(file.path(outpath, "library_depth.pdf"), depth_plots)
  }
  if (is.na(depth_filter)) {
    before_plot <- plot_library_depth(dds, "Library Depth Prefiltering", bins = 100)
    # ggsave(file.path(outpath, "library_depth.pdf"), before_plot)
  }
  
  # Take a look at library size
  message("Checking library size")
  if (!is.na(size_filter)) {
    before_plot <- plot_library_size(dds_before, "Library Size Prefiltering", bins = 10)
    keep <- colSums(as.data.frame(normal(dds))) >= size_filter
    dds <- dds[,keep]
    after_plot <- plot_library_size(dds, "Library Size Postfiltering", bins = 10)
    size_plots <- plot_grid(before_plot, after_plot, ncol = 2)
    # ggsave(file.path(outpath, "library_size.pdf"), size_plots)
  }
  if (is.na(size_filter)) {
    before_plot <- plot_library_size(dds, "Library Size Prefiltering", bins = 10)
    # ggsave(file.path(outpath, "library_size.pdf"), before_plot)
  }

  # Take a look at percent of genes detected
  message("Checking percent of genes detected")
  if (!is.na(percent_filter)) {
    before_plot <- plot_percent_genes_detected(dds_before, "Percent of Genes Detected Prefiltering")
    keep <- percentGenesDetected(dds) >= percent_filter
    dds <- dds[keep,]
    after_plot <- plot_percent_genes_detected(dds, "Percent of Genes Detected Postfiltering")
    percent_plots <- plot_grid(before_plot, after_plot, ncol = 2)
    # ggsave(file.path(outpath, "percent_genes_detected.pdf"), percent_plots)
  }
  if (is.na(percent_filter)) {
    before_plot <- plot_percent_genes_detected(dds, "Percent of Genes Detected Prefiltering")
    # ggsave(file.path(outpath, "/percent_genes_detected.pdf"), before_plot)
  }

  # save all the plots together
  all_plots <- plot_grid(depth_plots, size_plots, percent_plots, nrow = 3)
  ggsave(file.path(outpath, "preprocessing_plots.pdf"), all_plots)
  
  # Make a tree
  message("Checking hierarchical clustering")
  counts <- dds %>% normal() %>% t() %>% scale()
  sampleTree = hclust(dist(counts))
  ggtree(sampleTree) + 
    geom_tiplab(size = 2, hjust = -0.1) + 
    theme_tree2() +
    labs(title = "Sample clustering to detect outliers", 
         subtitle = "",
         x = "",
         y = "")
  ggsave(file.path(outpath, "sample_outliers.pdf"), dpi = 300)

  # PCA Plot
  message("Checking PCA")
  # get the variance stabilized counts
  counts <- dds %>% vst() %>% assay() %>% t() %>% scale()
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
  counts <- dds %>% assay() %>% t() %>% scale()
  sampleTree = hclust(dist(counts))
  ggtree(sampleTree) + 
    geom_tiplab(size = 2, hjust = -0.1) + 
    theme_tree2() +
    labs(title = "Sample clustering to detect outliers", 
         subtitle = "",
         x = "",
         y = "")
  ggsave(file.path(outpath, "sample_outliers.pdf"), dpi = 300)

  # PCA Plot
  message("Checking PCA")
  # get the variance stabilized counts
  counts <- dds %>% vst() %>% assay() %>% t() %>% scale()
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

  return(dds)
}

# Wilcoxan Ranked Sum testing on genes in two summarized experiments
# Arguments:
#   dds: DESeq2 object
#   condition: column of interest
#   controls: control groups
#   pCutoff: pvalue cutoff for volcano plot
#   FCcutoff: fold change cutoff for volcano plot
#   ...: additional arguments to pass to wilcox.test
# Outputs:
#   data frame of p-values and adjusted p-values
gene_wilcox_test <- function(
  dds, outpath,
  condition,
  pCutoff = 0.05, FCcutoff = 0.5,
  ...) {
  require(DESeq2)

  # make the outpath
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Extract count data from the DESeq object
  count_data_raw <- counts(dds, normalize = T) %>% 
    t() %>% 
    as.data.frame()
  count_data <- log2(count_data_raw + 1)
  counts_data_raw <- normalize_counts(dds, method = 'log2-mor')
  
  # Extract metadata from the DESeq object
  meta <- colData(dds) %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(condition))
  
  # get the conditions
  conditions <- unique(meta[, condition])
  
  # Run the wilcoxon test for each gene using map
  test_res <- map(
    count_data, 
    function(x) {
      # run the test
      test <- wilcox.test(
        x[meta[, condition] == conditions[1]],
        x[meta[, condition] == conditions[2]],
        ...
        )
      
      # get the pvalue
      pvalue <- test$p.value

      # get the basemean
      basemean <- 2^mean(x)

      # get the log2FC
      log2FoldChange <- mean(x[meta[, condition] == conditions[2]]) - mean(x[meta[, condition] == conditions[1]])
      
      # return the test
      return(list(basemean = basemean, log2FoldChange = log2FoldChange, pvalue = pvalue))
    }
    )
  names(test_res) <- colnames(count_data)


  # Extract p-values and adjust for multiple testing using the Benjamini-Hochberg method
  res <- list_of_lists_to_df(test_res)
  res <- res %>% mutate(padj = p.adjust(pvalue, method="BH"))

  # wrote the results to a csv
  write.csv(res, file.path(outpath, "wilcox_results.csv"))

  # make a volcano plot
  volcanoP <- EnhancedVolcano(
    res, lab=rownames(res), 
    x = 'log2FoldChange', y = 'pvalue', 
    title = 'Volcano Plot', subtitle = '',
    pCutoff = pCutoff, FCcutoff = FCcutoff
    )
  ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # get the fc list
  fc <- get_fc_list(res, 'log2FoldChange')

  # run enrichment
  gse <- rna_enrichment(fc, outpath)
  
  # View results
  return(res)
}

# Function to run limma voom on a summarized experiment object with TMM normalization
# Arguments:
#   dds: DESeq2 object
#   condition: vector of conditions
#   controls: vector of control groups
#   ...: additional arguments to pass to voom
# Outputs:
#   data frame of p-values and adjusted p-values
run_limma <- function(
    se, outpath,
    condition, controls = NULL, 
    pCutoff = 0.05, FCcutoff = 0.5,
    ...) {
    # make the directory if it doesn't exist
    dir.create(outpath, showWarnings = F, recursive = T)
    
    # Convert the counts matrix to a matrix if it's not already
    counts <- assay(se)
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }

    condition <- se[[condition]]

    if (!is.null(controls)) {
        controls <- se[[controls]]
    }

    # Normalize the counts matrix using TMM
    require(edgeR)
    counts <- DGEList(counts = counts, group = condition)
    counts <- calcNormFactors(counts, method = "TMM")
    norm_counts <- cpm(counts, log = TRUE)

    # Create the design matrix
    if (is.null(controls)) {
        design <- model.matrix(~ condition)
    } else {
        design <- model.matrix(~ condition + controls)
    }

    # Perform differential expression analysis
    fit <- lmFit(norm_counts, design)
    fit <- eBayes(fit)

    # Get the differential expression results
    results <- topTable(fit, coef = 2, number = Inf, ...)

    # make a volcano plot
    volcanoP <- EnhancedVolcano(
        results, lab=rownames(results), 
        x = 'logFC', y = 'P.Value', 
        title = 'Volcano Plot', subtitle = '',
        pCutoff = pCutoff, FCcutoff = FCcutoff
    )
    ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

    # save the results
    write.csv(results, file.path(outpath, "limma_results.csv"))

    # get the fc list
    fc <- get_fc_list(results, 'logFC')

    # run enrichment
    gse <- rna_enrichment(fc, outpath)

    return(results)
}

#  Calculate Correlations
#  Arguments:
#    dds: DESeq2 object
#    condition: column of interest
#    method: correlation method to use
#  Outputs:
#    Correlation results for each gene
calculate_correlations <- function(dds, condition, normalize = 'mor', method="pearson", ...) {
  # Extract the condition vector
  condition <- colData(dds)[, condition] %>% as.numeric()
  
  # Extract the gene expression data
  gene_expression <- normalize_counts(dds, method = normalize)
  
  # Apply the cor() function to each row of the gene_expression matrix
  gene_expression_correlations <- apply(gene_expression, 1, function(x) {
    correlation <- cor.test(x, condition, method=method, ...)
    return(list(correlation=correlation$estimate, pvalue=correlation$p.value))
  })
  
  # Convert the correlations to a data frame
  gene_expression_correlations_df <- data.frame(row.names = rownames(dds),
                                                correlation=sapply(gene_expression_correlations, function(x) x$correlation),
                                                pvalue=sapply(gene_expression_correlations, function(x) x$pvalue)
  )
  # Return the data frame of gene expression correlations
  return(gene_expression_correlations_df)
}

#  Run Simple Deseq analysis
# Arguments:
#   dds: DESeq2 object
#   outpath: path to push results to
#   contrast: contrast to run
#   ...: additional arguments to pass to DESeq2
# Outputs:
#   results of differential expression analysis
run_deseq <- function(
  # main inputs
  dds, outpath,
  # additional arguments for DESeq2 
  contrast = NA, ...,
  # additional arguments for plotting
  pvalue = "padj", pCutoff = 0.05, FCcutoff = 0
  ) {
  require(DESeq2)
  require(SummarizedExperiment)
  require(tidyverse)
  require(ggplot2)

  # This function will run differential expression on a premade dds object.
  # Gives most metrics you'll need, and also returns the results
  dir.create(outpath, showWarnings = F, recursive = T)
  message("Running DESeq2")
  dds <- DESeq(dds, ...)

  # Get results names
  res <- results(dds)
  
  if (is.vector(contrast)) {
    name <- paste0(contrast[1], '__', contrast[2], '_vs_', contrast[3])
    res <- results(dds, contrast=contrast)
    resLFC <- lfcShrink(dds, coef=length(resultsNames(dds)), type="apeglm")
    pdf(file.path(outpath, "MAplot.pdf"))
    DESeq2::plotMA(resLFC)
    dev.off()

    # make volcano plot
    # volcanoP <- EnhancedVolcano(
    #   res, lab=rownames(res), 
    #   x = 'log2FoldChange', y = 'pvalue', 
    #   title = name, subtitle = '', 
    #   pCutoff = pCutoff, FCcutoff = FCcutoff
    #   )
    volcanoP <- plot_volcano(
      res, 
      title = name,
      color = pvalue,
      pCutoff = pCutoff
    )
    ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

    # make gene list
    geneList <- get_fc_list(res)
    # gse <- rna_enrichment(geneList, outpath)
    gse_list <- gsea_analysis(geneList, outpath)


    message("Results Summary:")
    summary(res)
    
    message("Writing results to outpath")
    write.csv(res, file.path(outpath, "deseq_results.csv"))
    saveRDS(dds, file = file.path(outpath,"dds.rds"))
    return(res)
  }

  name <- resultsNames(dds)[length(resultsNames(dds))]
  res <- results(dds, contrast=contrast)
  resLFC <- lfcShrink(dds, coef=name, type="apeglm")
  pdf(file.path(outpath, "MAplot.pdf"))
  DESeq2::plotMA(resLFC)
  dev.off()

  # make volcano plot
  volcanoP <- EnhancedVolcano(
    res, lab=rownames(res), 
    x = 'log2FoldChange', y = pvalue, 
    title = name, subtitle = '',
    pCutoff = pCutoff, FCcutoff = FCcutoff
    )
  ggsave(file.path(outpath, "volcanoPlot.pdf"), volcanoP)

  # make a heatmap
  sign_genes <- rownames(res)[res$pvalue < pCutoff & abs(res$log2FoldChange) > FCcutoff]
  heatmapP <- plot_gene_heatmap(
    dds[], 
    title = name,
    annotations = contrast[1],
    normalize = "vst"
    )
  pdf(file.path(outpath, "dge_heatmap.pdf"))
  print(heatmapP)
  dev.off()

  # make gene list
  geneList <- get_fc_list(res)
  # gse <- rna_enrichment(geneList, outpath)
  gse_list <- gsea_analysis(geneList, outpath)
  
  message("Results Summary:")
  summary(res)
  
  message("Writing results to outpath")
  write.csv(res, file.path(outpath, "deseq_results.csv"))
  saveRDS(dds, file = file.path(outpath,"dds.rds"))
  return(res)
}

# OVR Deseq Results Function
# Arguments:
#   dds: DESeq2 object
#   column: column of interest for OVR analysis
#   outpath: path to push results to
# Outputs:
#   OVR results for each level of column of interest
ovr_deseq_results <- function(dds, column, outpath, ...) {
  require(DESeq2)
  require(clusterProfiler)
  require(DOSE)
  require(enrichplot)

  dir.create(outpath, showWarnings = F, recursive = T)

  # This function takes in a DDS or SE object, the condition column of interest
  # and the out directory path to push results to. It will run OVR differential
  # expression analysis on each level within the condition column.
  counts <- assay(dds)
  lvls <- levels(colData(dds)[,column])
  cond <- as.character(colData(dds)[,column])

  # loop over condition levels in a one versus rest manner
  # ideally this for loop could be an apply statement with a custom function, 
  # but it works for now so oh well lol
  list_out <- list()
  for (lvl in lvls) {
    print(paste0('Testing ', column, ' ', lvl, ' versus rest'))
    # Set our OVR analysis
    cond_ <- cond
    cond_[cond_ != lvl] <- 'rest'
    
    # create temporary dds objects for analysis
    dds_ <- DESeqDataSetFromMatrix(counts, colData <- DataFrame(condition = as.factor(cond_)), 
                                   design <- ~ condition)
    dds_ <- DESeq(dds_)
    res <- results(dds_, contrast = c('condition', lvl, 'rest'))

    path <- file.path(outpath, paste0(column, '__',lvl,'_v_rest'))
    dir.create(path, showWarnings = F, recursive = T)

    saveRDS(dds_, file=paste0(path, '/deseqDataset_', column,'__',lvl,'_v_rest.rdata'))
    write.csv(res, file=paste0(path, '/dge_results_', column,'__',lvl,'_v_rest.csv'))
    
    volcanoP <- EnhancedVolcano(res, lab=rownames(res), x = 'log2FoldChange', y = 'pvalue', title = paste0(column,'__',lvl,'_v_rest'), subtitle = '', ...)
    ggsave(paste0(path, '/volcanoPlot_', column,'__',lvl,'_v_rest.pdf'), volcanoP)

    geneList <- get_fc_list(res)
    gse <- gsea_analysis(geneList, path)



    list_out[[lvl]] <- res
    }
  return(list_out)
}

# Function to remove NA variables from a summarized experiment object
# Arguments:
#   se: SummarizedExperiment object
#   columns: list of columns to remove NAs from
# Outputs:
#   DESeqDataSet object with NAs removed
remove_na_variables <- function(se, columns) {
  # Get the colData from the summarized experiment object
  col_data <- colData(se)
  assay_data <- assay(se)
  
  # Find the variables with missing values
  na_samples <- rownames(col_data)[rowSums(is.na(col_data[, columns])) > 0]

  # Remove the variables with missing values from the colData
  col_data <- col_data[!(rownames(col_data) %in% na_samples),]
  
  # Update the DESeqDataSet object with the modified colData
  new_se <- make_se(assay_data, col_data)

  # Return the new DESeqDataSet object
  return(new_se)
}

# Function to take in dds object and deseq results and return summary of results
# Arguments:
#   results: results of differential expression analysis
#   padj_cutoffs: list of padj cutoffs to use
#   pvalue_cutoffs: list of pvalue cutoffs to use
# Outputs:
#   summary of results
summarize_deseq_experiment <- function(
  results, 
  padj_cutoffs = c(0.05, 0.1, 0.2), pvalue_cutoffs = c(0.01, 0.05, 0.1), 
  logFC_cutoff = 0
  ) {
  # summary of the results at different padj cutoffs
  padj_summary <- data.frame()
  for (padj_cutoff in padj_cutoffs) {
    padj_summary <- rbind(
      padj_summary, data.frame(sign_cutoff = paste0('padj < ', padj_cutoff), 
      n_genes = filter(as.data.frame(results), padj < padj_cutoff) %>% nrow(), 
      n_up = filter(as.data.frame(results), padj < padj_cutoff & log2FoldChange > logFC_cutoff) %>% nrow(),
      n_down = filter(as.data.frame(results), padj < padj_cutoff & log2FoldChange < logFC_cutoff) %>% nrow())
      )
  }

  # summary of the results at different pvalue cutoffs
  pvalue_summary <- data.frame()
  for (pvalue_cutoff in pvalue_cutoffs) {
    pvalue_summary <- rbind(
      pvalue_summary, 
      data.frame(sign_cutoff = paste0('pvalue < ', pvalue_cutoff), n_genes = filter(as.data.frame(results), pvalue < pvalue_cutoff) %>% nrow(),
      n_up = filter(as.data.frame(results), pvalue < pvalue_cutoff & log2FoldChange > logFC_cutoff) %>% nrow(),
      n_down = filter(as.data.frame(results), pvalue < pvalue_cutoff & log2FoldChange < -logFC_cutoff) %>% nrow())
      )
  }

  # combine the summaries
  summary <- rbind(padj_summary, pvalue_summary)
  return(summary)
}

# Function to run deseq on a variety of conditions
# Arguments:
#   dds: DESeq2 object
#   conditions: list of conditions to run deseq on
#   controls: list of controls to run deseq on
#   outpath: path to push results to
#   ...: additional arguments to pass to deseq for volcano plots
# Outputs:
#   results of differential expression analysis
deseq_analysis <- function(dds, conditions, controls, outpath, ...) {
  # make directory
  dir.create(outpath, showWarnings = F, recursive = T)

  # make list to store results
  analysis_list <- list()
  
  # make summary dataframe 
  summary_df <- data.frame()
  
  # make design matrix
  for (condition in conditions) {
    # make directory for condition
    dir.create(file.path(outpath, condition), showWarnings = F, recursive = T)

    # make new dds object with no NAs
    dds_ <- remove_na_variables(dds, c(controls, condition))

    # make a stats table of the conditions and controls
    df_stats <- as.data.frame(colData(dds_))
    stats_table <- stats_table(df_stats, condition, vars = controls)
    # save the stats table
    write.csv(stats_table, file.path(outpath, condition, paste0(condition, '_stats_table.csv')))

    # check if condition is a factor
    if (!is.factor(colData(dds_)[, condition])) {
      colData(dds_)[, condition] <- as.factor(colData(dds_)[, condition])
      print(paste0('Converting ', condition, ' to factor'))
    }

    input_ <- paste0(append(controls, condition), collapse = ' + ')
    design_matr <- as.formula(paste0("~ ", input_))
    dds_ <- DESeqDataSet(dds_, design = design_matr)
    levels <- levels(colData(dds_)[, condition])

    message(paste0('Running Analysis on ', condition))
    sink(file.path(outpath, condition, paste0(condition, '_deseq_analysis_summary.log')))
    # if there are more than 2 levels, run OVR analysis
    if (length(levels) > 2) {
      res <- ovr_deseq_results(dds_, condition, file.path(outpath, condition), ...)

      # # summarize res
      summary <- lapply(res, summarize_deseq_experiment)

      # make condition for the res
      # condition <- names(res)
      for (i in 1:length(summary)) {
        summary[[i]]$condition <- names(res)[i]
      }
      # summary <- lapply(summary, function(x) {x$condition <- condition; return(x)})
      summary <- do.call(rbind, summary)

      # reorder columns so condition is first
      summary <- summary[, c(ncol(summary), 1:(ncol(summary)-1))]

      # add summary to dataframe
      summary_df <- rbind(summary_df, summary)

      # add results to list
      analysis_list[[condition]] = res
    }

    # if there are only 2 levels, run normal deseq
    if (length(levels) == 2) {
      # flip levels if they are 0 and 1
      # if(setequal(levels, c('0', '1'))) {levels <- c('1', '0')}

      # get contrast and run deseq
      contrast <- c(condition, levels[2], levels[1])
      res <- run_deseq(dds_, file.path(outpath, condition), contrast = contrast, ...)

      # summarize res
      summary <- summarize_deseq_experiment(res, padj_cutoffs = c(0.05, 0.1, 0.2), pvalue_cutoffs = c(0.01, 0.05, 0.1))
      summary$condition <- condition

      # reorder columns so condition is first
      summary <- summary[, c(ncol(summary), 1:(ncol(summary)-1))]

      # add summary to dataframe
      summary_df <- rbind(summary_df, summary)

      # add results to list
      analysis_list[[condition]] = res
    }

    # if there is only one level, skip
    if (length(levels) == 1) {
      message(paste0('Skipping ', condition, ' because there is only one level'))
    }

    # if there are no levels, skip
    if (length(levels) == 0) {
      message(paste0('Skipping ', condition, ' because there are no levels'))
    }

  sink()
  }
  # save summary dataframe
  write.csv(summary_df, file.path(outpath, 'deseq_analysis_summary.csv'), row.names = FALSE)
  return(analysis_list)
}

# Function to run statistical tests on a DESeq2 object using rstatix functions
# Arguments:
#   formula: formula to use for the test
#   dds: DESeq2 object
#   rstatix_test: rstatix test to use
#   ...: additional arguments to pass to rstatix_test
# Outputs:
#   data frame of p-values and adjusted p-values
test_dds <- function(formula, dds, rstatix_test = wilcox_test, ...) {
    require(DESeq2)
    message("Extracting data from DESeq object")
    # Extract count data from the DESeq object
    count_data_raw <- assay(dds) %>% 
      t() %>% 
      as.data.frame()
    count_data <- log2(count_data_raw + 1)
    
    message("Extracting metadata from DESeq object")
    # Extract metadata from the DESeq object
    meta <- colData(dds) %>% 
      as.data.frame() %>%
      dplyr::select(all_of(attr(terms(formula), "term.labels")))

    message("Running wilcoxon test")
    # run the test for all the genes using map
    test_res <- map(
      count_data,
      function(x) {
        # combine the metadata and count data on the gene
        data <- meta %>% mutate(x = x)

        # add the gene to the formula as the response using reformulate
        formula <- reformulate(
          response = 'x',
          termlabels = attr(terms(formula), "term.labels")
          )

        # get the condition as the last term in the formula
        condition <- attr(terms(formula), "term.labels")[length(attr(terms(formula), "term.labels"))]
        conditions <- unique(meta[, condition])
        
        # run the test
        test <- rstatix_test(
          data,
          formula,
          ...
          )

        # get the basemean
        basemean <- 2^mean(x)

        # get the log2FC
        log2FC <- mean(x[meta[, condition] == conditions[2]]) - mean(x[meta[, condition] == conditions[1]])

        # add the basemean and log2FC to the test
        test$basemean <- basemean
        test$log2FC <- log2FC
        
        # return the test
        return(test)
      }
      )

      # add the gene names
      names(test_res) <- colnames(count_data)

    # convert to dataframe
    df_out <- list_of_lists_to_df(test_res) %>%
      mutate(
        pvalue = p,
        padj = p.adjust(p, method = "fdr")
      ) %>%
      dplyr::select(
        basemean,
        log2FC,
        pvalue,
        padj
      )

    # return the dataframe
    return(df_out)
}






















#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
