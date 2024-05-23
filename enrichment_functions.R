###########################################################################
#
#                            enrichment_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: enrichment_functions

#======================== LIBRARIES ========================#
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(enrichplot))
suppressMessages(library(forcats))
suppressMessages(library(ggplot2))
suppressMessages(library(msigdbr))
suppressMessages(library(ggpubr))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/converting_functions.R")


#======================== CODE ========================#

# function to get the fc list
# Arguments:
#   res: data frame of results
#   fc_col: column of fc values
#   names: column of names (if null default to rownames)
# Outputs:
#   vector of sorted fc values
get_fc_list <- function(res, fc_col = "log2FoldChange", names = NULL) {
    if (is.null(names)) {
        res[, "rownames"] <- rownames(res)
        names <- "rownames"
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
  
  image_type = "pdf", ontology = "ALL", 
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

  gse <- do.call(enricher_function, list(geneList, org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf, ...))
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)
  saveRDS(gse, file.path(outpath, "enrichment_results.rds"))
  save_gse(gse, outpath, image_type = image_type)
  return(gse)
}

# Function to save and plot gse object
# Arguments:
#   gse: gse object
#   outpath: path to save to
#   ...: additional arguments to pass to ggsave
# Outputs:
#   saves the gse object to the outpath
save_gse <- function(gse, outpath, ...) {
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  require(enrichplot)
  require(ggplot2)

  # save the gse object
  write.csv(gse@result, file.path(outpath, "enrichment_results.csv"), quote = TRUE, row.names = FALSE)
  write.csv(filter(gse@result, qvalue < 0.1), file.path(outpath, "enrichment_results_sig.csv"), quote = TRUE, row.names = FALSE)
  saveRDS(gse, file.path(outpath, "enrichment_results.rds"))

  tryCatch({
  # Dotplot
  gseDot <- enrichplot::dotplot(gse, showCategory = 20) + ggtitle("Enrichment Dotplot")
  ggsave(file.path(outpath, paste0("dotplot.pdf")), gseDot, ...)
  }, error = function(e) {
    warning("Dotplot GSEA Plots Failed")
  })

  tryCatch({
    # cnetplot
    cnet <- cnetplot(gse, node_label="category", cex_label_gene = 0.8)
    ggsave(file.path(outpath, paste0("cnetplot.pdf")), cnet, ...)
  }, error = function(e) {
    warning("Cnetplot GSEA Plots Failed")
  })

  tryCatch({
    # Barplot data
    gse_bar <- gse %>%
      as.data.frame() %>%
      group_by(sign(NES)) %>%
      arrange(qvalue) %>%
      slice(1:10)

    # Barplot
    gseBar <- ggplot(gse_bar, aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
      geom_col(orientation = "y") +
      scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) +
      labs(title="Enrichment Barplot", y = NULL) +
      theme_classic2()
    ggsave(file.path(outpath, paste0("barplot.pdf")), gseBar, ...)
  }, error = function(e) {
    warning("GSEA Barplot Failed")
  })

  tryCatch({
    # Platelet Termed Barplot
    p_terms <- plot_enrichment_terms(gse, terms2plot = c("inflam", "immune", "plat", "coag"))
    ggsave(file.path(outpath, paste0("barplot_terms.pdf")), p_terms, ...)
  }, error = function(e) {
    warning("GSEA Term Specific Barplot Failed")
  })

  tryCatch({
    p_ridge <- ridgeplot(gse, showCategory = 15)
    ggsave(file.path(outpath, paste0("ridgeplot.pdf")), p_ridge, ...)
  }, error = function(e) {
    warning("RidgePlot GSEA Failed")
  })

    tryCatch({
    p_heat <- heatplot(gse, showCategory = 10)
    ggsave(file.path(outpath, paste0("heatplot.pdf")), p_heat, ...)
  }, error = function(e) {
    warning("Heatplot GSEA Failed")
  })

}

#  Run a gsea analysis
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

  msigdb_category = "H", # I"m going to depricate this
  ontology = "ALL"
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
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # Run GSEA on a few genesets
  msigdb <- msigdbr(species = "Homo sapiens")
  gse_go <- gseGO(geneList, org.Hs.eg.db, keyType = keyType, ont = ontology, pvalueCutoff = Inf)

  H_t2g <- msigdb %>%
    filter(gs_cat == "H") %>%
    dplyr::select(gs_name, gene_symbol)
  gse_h <- GSEA(geneList, TERM2GENE = H_t2g, pvalueCutoff = Inf)

  reactome_t2g <- msigdb %>%
    filter(gs_cat == "C2" & gs_subcat == "CP:REACTOME") %>%
    dplyr::select(gs_name, gene_symbol)
  gse_reactome <- GSEA(geneList, TERM2GENE = reactome_t2g, pvalueCutoff = Inf)

  # gse list to loop over
  gse_list <- list(
    GO = gse_go,
    H = gse_h,
    REACTOME = gse_reactome
    )

  for (idx in seq_along(gse_list)) {
    # get the gse and the name
    gse <- gse_list[[idx]]
    name <- names(gse_list)[idx]
    save_gse(gse, file.path(outpath, name))
  }
  return(gse_list)
}