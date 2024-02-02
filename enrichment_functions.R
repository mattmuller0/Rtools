###########################################################################
#
#                            enrichment_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: enrichment_functions

#======================== LIBRARIES ========================#
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(forcats)
library(ggplot2)
library(msigdbr)
library(ggpubr)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')


#======================== CODE ========================#

# function to get the fc list
# Arguments:
#   res: data frame of results
#   fc_col: column of fc values
#   names: column of names (if null default to rownames)
# Outputs:
#   vector of sorted fc values
get_fc_list <- function(res, fc_col = 'log2FoldChange', names = NULL) {
    if (is.null(names)) {
        res[, 'rownames'] <- rownames(res)
        names <- 'rownames'
    }

    # get the fc values
    fc <- res %>% 
        as.data.frame() %>%
        # get the fc values
        dplyr::select(!!sym(fc_col)) %>% 
        # make a named vector
        unlist() %>%
        # name the vector
        set_names(res[, names]) %>%
        # order the vector
        sort(decreasing = T)
    
    # return the fc
    return(fc)
}

#  Run simple enrichment with enrichGO or gseGO
#  Arguments:
#    geneList: list of genes to run enrichment on
#    outpath: path to push results to
#    keyType: key type for gene list
#    enricher: enrichment method to use (default is gseGO)
#    image_type: type of image to save (default is pdf)
#    ...: additional arguments to pass to enricher
#  Outputs:
#    Enrichment results for each level of column of interest
rna_enrichment <- function(
  geneList, outpath, 
  keyType = NULL, enricher_function = NULL, 
  
  image_type = 'pdf', ontology = "ALL", 
  terms2plot = c("inflam", "immune", "plat"), 
  ...
  ) {
  require(SummarizedExperiment)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  require(enrichplot)

  # get the type of the gene list if not specified
  if (is.null(keyType)) {keyType <- detect_gene_id_type(names(geneList), strip = TRUE)}

  # Output directory
  dir.create(outpath, showWarnings = F, recursive = T)

  # Get enrichment type and go on
  if (is.null(enricher_function)) {
    message("No enrichment function specified, defaulting to gseGO")
    enricher_function <- gseGO
    }

  # gse <- enricher_function(geneList, org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf, ...)
  gse <- do.call(enricher_function, list(geneList, org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf, ...))
  
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)


  tryCatch(
    {
    # Dotplot
    gseDot <- enrichplot::dotplot(gse, showCategory = 20) + ggtitle("Enrichment Dotplot");
    ggsave(file.path(outpath, paste0('gseDotPlot.pdf')), gseDot)
    }, 

    # error handling
    error = function(e) {warning("Failed to make dotplot")}
  )

  require(forcats)
  require(ggplot2)
  require(enrichplot)

  tryCatch(
    {
    # Barplot data
    gse_bar <- gse %>%
      group_by(sign(NES)) %>%
      arrange(qvalue) %>%
      slice(1:10)

    # Barplot
    gseBar <- ggplot(gse_bar, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=20) + 
      geom_col(orientation='y') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
      labs(title='Enrichment Barplot') +
      theme_classic2() + ylab(NULL)
    ggsave(file.path(outpath, paste0('gseBarPlot.', image_type)), gseBar)

    # term specific barplot
    p_terms <- plot_enrichment_terms(gse, terms2plot = terms2plot)
    ggsave(file.path(outpath, paste0('gseBarPlot_terms.', image_type)), p_terms)
    },

    # error handling
    error = function(e) {warning("Failed to make barplot")}
  )
  return(gse)
}

#  Run simple enrichment with enrichGO or gseGO
#  Arguments:
#    geneList: list of genes to run enrichment on
#    outpath: path to push results to
#    keyType: key type for gene list
#    enricher: enrichment method to use (default is gseGO)
#    image_type: type of image to save (default is pdf)
#    ...: additional arguments to pass to enricher
#  Outputs:
#    Enrichment results for each level of column of interest
gsea_analysis <- function(
  geneList, outpath, 
  keyType = NULL,
  terms2plot = c("inflam", "immune", "plat"),
  genes2plot = c('SREB', 'VWF', 'CXLC5', 'GCA'),

  msigdb_category = "H", ontology = "ALL"
  ) {
  require(SummarizedExperiment)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  require(enrichplot)
  require(forcats)
  require(ggplot2)
  require(msigdbr)

  # get the type of the gene list if not specified
  if (is.null(keyType)) {keyType <- detect_gene_id_type(names(geneList), strip = TRUE)}

  # Output directory
  dir.create(outpath, showWarnings = F, recursive = T)

  # prep the gene list for kegg
  entrez_names <- map_gene_ids(names(geneList), from = keyType, to = "ENTREZID", remove_missing = TRUE)
  entrez_gL <- geneList
  names(entrez_gL) <- entrez_names

  # get the hallmark geneset
  msigdbr_H <- msigdbr(species = "Homo sapiens", category = msigdb_category)
  H_t2g <- msigdbr_H %>% 
    dplyr::select(gs_name, entrez_gene)

  # Run GSEA on a few genesets
  gse_go <- gseGO(geneList, org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf)
  # gse_kegg <- gseKEGG(entrez_gL, organism = 'hsa', pvalueCutoff = Inf)
  # gse_msigdb <- GSEA(entrez_gL, TERM2GENE = H_t2g, minGSSize = 10, maxGSSize = 500, pvalueCutoff = Inf)

  # gse list to loop over
  gse_list <- list(
    GO = gse_go,
    KEGG = gse_kegg,
    msigdb = gse_msigdb
  )

  for (idx in 1:length(gse_list)) {
    # get the gse and the name
    gse <- gse_list[[idx]]
    name <- names(gse_list)[idx]
    
    write.csv(gse@result, file.path(outpath, paste0(name, "_enrichment_results.csv")), quote = TRUE, row.names = FALSE)

    # dotplot with error catching
    tryCatch(
      {
      # Dotplot
      gseDot <- enrichplot::dotplot(gse, showCategory = 20) + ggtitle("Enrichment Dotplot");
      ggsave(file.path(outpath, paste0(name, '_DotPlot.pdf')), gseDot)
      }, 

      # error handling
      error = function(e) {warning("Failed to make dotplot")}
    )
    
    tryCatch(
      {
        gse_bar <- gse %>% 
          as.data.frame() %>%
          group_by(sign(NES)) %>%
          arrange(qvalue) %>%
          slice(1:10)

        # Barplot
        gseBar <- ggplot(gse_bar, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=20) + 
          geom_col(orientation='y') + 
          scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
          labs(title='Enrichment Barplot') +
          theme_classic2() + ylab(NULL)
        
        ggsave(file.path(outpath, paste0(name, '_BarPlot.pdf')), gseBar)

      },

      # error handling
      error = function(e) {warning("Failed to make barplot")}
  )

  # terms barplot with error catching
    tryCatch(
      {
      # term specific barplot
      p_terms <- plot_enrichment_terms(gse, terms2plot = terms2plot)
      ggsave(file.path(outpath, paste0(name, '_BarPlot_terms.pdf')), p_terms)
      },

      # error handling
      error = function(e) {warning("Failed to make terms barplot")}
    )

    # genes barplot with error catching
    tryCatch(
      {
      # term specific barplot
      p_genes <- plot_enrichment_terms(gse, genes2plot = genes2plot)
      ggsave(file.path(outpath, paste0(name, '_BarPlot_genes.pdf')), p_genes)
      },

      # error handling
      error = function(e) {warning("Failed to make genes barplot")}
    )
  }
  return(gse_list)
}












