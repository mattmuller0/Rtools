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
# suppressMessages(library(AnnotationDbi))
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

#' function to get the fc list
#' Arguments:
#'   res: data frame of results
#'   fc_col: column of fc values
#'   names: column of names (if null default to rownames)
#' Outputs:
#'   vector of sorted fc values
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

#'  Run simple enrichment with enrichGO or gseGO
#'  Arguments:
#'    geneList: list of genes to run enrichment on
#'    outpath: path to push results to
#'    keyType: key type for gene list
#'    enricher: enrichment method to use (default is gseGO)
#'    image_type: type of image to save (default is pdf)
#'    ...: additional arguments to pass to enricher
#'  Outputs:
#'    Enrichment results for each level of column of interest
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

#' Function to save and plot gse object
#' Arguments:
#'   gse: gse object
#'   outpath: path to save to
#'   ...: additional arguments to pass to ggsave
#' Outputs:
#'   saves the gse object to the outpath
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
      slice(1:10) %>%
      mutate(
        Description = gsub("^(REACTOME_|GO_|HALLMARK_)", "", Description),
        Description = gsub("_", " ", Description),
        Description = factor(stringr::str_wrap(Description, 40))
        )

    # Barplot
    gseBar <- ggplot(gse_bar, aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
      geom_col(orientation = "y") +
      scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) +
      labs(title="Enrichment Barplot", y = NULL) +
      theme_classic2()
    ggsave(file.path(outpath, paste0("barplot_all.pdf")), gseBar, ...)
  }, error = function(e) {
    warning("GSEA Barplot Failed")
  })

  tryCatch({
    # Barplot data
    gse_bar_bp <- gse %>%
      as.data.frame() %>%
      filter(ONTOLOGY == "BP") %>%
      group_by(sign(NES)) %>%
      arrange(qvalue) %>%
      slice(1:10) %>%
      mutate(
        Description = gsub("^(REACTOME_|GO_|HALLMARK_)", "", Description),
        Description = gsub("_", " ", Description),
        Description = factor(stringr::str_wrap(Description, 40))
        )

    gseBar_bp <- ggplot(gse_bar_bp, aes(NES , fct_reorder(Description, NES), fill=qvalue)) +
      geom_col(orientation = "y") +
      scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) +
      labs(title="Enrichment Barplot", y = NULL) +
      theme_classic2()
    ggsave(file.path(outpath, paste0("barplot_BP.pdf")), gseBar_bp, ...)
  }, error = function(e) {
    warning("GSEA Barplot BP Failed")
  })

  tryCatch({
    # cnetplot
    cnet <- cnetplot(gse, node_label="category", cex_label_gene = 0.8)
    ggsave(file.path(outpath, paste0("cnetplot.pdf")), cnet, ...)
  }, error = function(e) {
    warning("Cnetplot GSEA Plots Failed")
  })

  tryCatch({
    # Platelet Termed Barplot
    p_terms <- plot_enrichment_terms(gse, terms2plot = c("inflam", "plat", "coag"), max_terms = 12)
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

#'  Run a gsea analysis
#'  Arguments:
#'    geneList: list of genes to run enrichment on
#'    outpath: path to push results to
#'    keyType: key type for gene list
#'    enricher: enrichment method to use (default is gseGO)
#'    image_type: type of image to save (default is pdf)
#'    ...: additional arguments to pass to enricher
#'  Outputs:
#'    Enrichment results for each level of column of interest
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

#' Function to do up and down overrepresentation analysis
#' Arguments
#'  - gene_dataframe [feature, direction] (output of getGenes)
#'  - method
#'  - padj_cutoff
#'  - max_pathways
stratified_ora <- function(
  gene_dataframe,
  outpath,
  method = "enrichGO",
  padj_cutoff = 0.05,
  max_pathways = 5,
  ...
  ) {
  require(clusterProfiler)
  require(enrichR)
  require(org.Hs.eg.db)

  # get the up and down genes
  up_genes <- gene_dataframe %>% filter(direction == "up") %>% pull(features)
  down_genes <- gene_dataframe %>% filter(direction == "down") %>% pull(features)

  # make sure the method is supported
  methods <- c("enrichGO", "groupGO")
  if (!(method %in% methods)) {
    stop("Method not supported. Please use one of: ", paste(methods, collapse = ", "))
  }

  enr_fn <- switch(method,
    "enrichGO" = function(x) enrichGO(x, org.Hs.eg.db, pvalueCutoff = Inf, ...),
    "groupGO" = function(x) groupGO(x, org.Hs.eg.db, pvalueCutoff = Inf, ...),
  )
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  out <- purrr::map_dfr(
    c("up", "down"), ~{
      if (.x == "up") {
        enr <- enr_fn(up_genes)
      } else {
        enr <- enr_fn(down_genes)
      }

      if (is.null(enr)) {
        message("No results found for ", .x)
        return(NULL)
      }
      save_gse(enr, file.path(outpath, .x))

      enr@result %>% 
        filter(p.adjust < padj_cutoff) %>%
        arrange(p.adjust) %>%
        slice(1:max_pathways) %>%
        mutate(direction = .x)
    }
  )
  write.csv(out, file.path(outpath, "ora_results.csv"), quote = TRUE, row.names = FALSE)
  out$signed <- ifelse(out$direction == "up", -log10(out$p.adjust), log10(out$p.adjust))
  p <- ggplot(out, aes(x = signed, y = fct_reorder(stringr::str_wrap(Description, 40), signed), fill = direction)) +
    geom_col() +
    labs(title = "ORA Results", x = "Signed -log10(padj)", y = NULL) +
    theme_classic2()
  ggplot2::ggsave(file.path(outpath, "ora_results.pdf"), p)

  return(out)
  }


#' Function to do up and down overrepresentation analysis
#' Arguments
#'  - gene_dataframe [feature, direction] (output of getGenes)
#'  - method
#'  - padj_cutoff
#'  - max_pathways
stratified_enrichr <- function(
  gene_dataframe,
  outpath,
  dbs = c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "WikiPathway_2023_Human", "GWAS_Catalog_2023", "Reactome_2022", "MSigDB Hallmark 2020"),
  padj_cutoff = 0.05,
  max_pathways = 5,
  ...
  ) {
  require(clusterProfiler)
  require(enrichR)
  require(org.Hs.eg.db)

  # get the up and down genes
  up_genes <- gene_dataframe %>% filter(direction == "up") %>% pull(features)
  down_genes <- gene_dataframe %>% filter(direction == "down") %>% pull(features)
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  out_dbs <- purrr::map(dbs, ~{
    dir.create(file.path(outpath, .x), showWarnings = FALSE, recursive = TRUE)
    enr_up <- enrichr(up_genes, .x)[[.x]] %>% mutate(direction = "up")
    enr_down <- enrichr(down_genes, .x)[[.x]] %>% mutate(direction = "down")
    enr <- bind_rows(enr_up, enr_down)
    if (nrow(enr) == 0) {
      message("No results found for ", .x)
      return(NULL)
    }
    write.csv(enr, file.path(outpath, .x, "enrichr_results.csv"), quote = TRUE, row.names = FALSE)

    sign_enr <- enr %>% 
      group_by(direction) %>%
      filter(Adjusted.P.value < padj_cutoff) %>%
      arrange(Adjusted.P.value) %>%
      slice(1:max_pathways) %>%
      mutate(signed = ifelse(direction == "up", -log10(Adjusted.P.value), log10(Adjusted.P.value)))
    write.csv(sign_enr, file.path(outpath, .x, "enrichr_results_sig.csv"), quote = TRUE, row.names = FALSE)

    p <- ggplot(sign_enr, aes(x = signed, y = fct_reorder(stringr::str_wrap(Term, 40), signed), fill = direction)) +
      geom_col() +
      labs(title = "ORA Results", x = "Signed -log10(padj)", y = NULL) +
      theme_classic2()
    ggplot2::ggsave(file.path(outpath, .x, "enrichr_results.pdf"), p)
    sign_enr
  })
  names(out_dbs) <- dbs

  return(out_dbs)
  }