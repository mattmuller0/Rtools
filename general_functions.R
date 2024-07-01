# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2023-01-11
# 
# Script Name: General Functions
# 
# Notes:
# Lots of fun functions I've made to work with.


# LOAD LIBRARIES ------------------------------------------
suppressMessages(library(tidyverse))

# CODE BLOCK ----------------------------------------------
#' Function to summarize results more generally
#' Arguments:
#'   results: results of differential expression analysis
#'   logFC_column: column to use for logFC
#'   pvalue_column: column to use for pvalue
#'   padj_column: column to use for padj
#'   padj_cutoffs: list of padj cutoffs to use
#'   pvalue_cutoffs: list of pvalue cutoffs to use
#'   logFC_cutoff: logFC cutoff to use
#' Outputs:
#'   summary of results
summarize_experiment <- function(
  results, 
  logFC_column = "log2FoldChange",
  pvalue_column = "pvalue",
  padj_column = "padj",
  pvalue_cutoffs = c(0.01, 0.05, 0.1), 
  padj_cutoffs = c(0.05, 0.1, 0.2), 
  logFC_cutoff = 0
  ) {
  # summary of the results at different padj cutoffs
  summary <- data.frame(
    variable = character(),
    p_cutoff = numeric(),
    fc_cutoff = numeric(),
    n_sig = numeric(),
    n_up = numeric(),
    n_down = numeric()
  )
  lapply(pvalue_cutoffs, function(pvalue_cutoff) {
    summary <<- summary %>% 
      add_row(
        variable = "pvalue",
        cutoff = pvalue_cutoff,
        fc_cutoff = logFC_cutoff,
        n_sig = sum(results[[pvalue_column]] < pvalue_cutoff),
        n_up = sum(results[[pvalue_column]] < pvalue_cutoff & results[[logFC_column]] > logFC_cutoff),
        n_down = sum(results[[pvalue_column]] < pvalue_cutoff & results[[logFC_column]] < -logFC_cutoff)
      )
  })
  lapply(padj_cutoffs, function(padj_cutoff) {
    summary <<- summary %>% 
      add_row(
        variable = "padj",
        cutoff = padj_cutoff,
        fc_cutoff = logFC_cutoff,
        n_sig = sum(results[[padj_column]] < padj_cutoff),
        n_up = sum(results[[padj_column]] < padj_cutoff & results[[logFC_column]] > logFC_cutoff),
        n_down = sum(results[[padj_column]] < padj_cutoff & results[[logFC_column]] < -logFC_cutoff)
      )
  })
  return(summary)
}

#' Function to get genes from a results table
#' Arguments:
#'  - results: results table
#'  - pval: p-value cutoff
#'  - log2fc: log2 fold change cutoff
#'  - gene_col: gene column name
#'  - pval_col: p-value column name
#'  - log2fc_col: log2 fold change column name
#' Returns:
#' - dataframe of genes separated by up and down
getGenes <- function(results, pval = 0.05, metric = 0, name_col = "rownames", pval_col = "padj", metric_col = "log2FoldChange") {
    if (name_col == "rownames") {results <- rownames_to_column(results, var = name_col)}
    up <- results %>%
        dplyr::filter(!!sym(pval_col) < pval & !!sym(metric_col) > metric) %>%
        dplyr::pull(!!sym(name_col))
    down <- results %>%
        dplyr::filter(!!sym(pval_col) < pval & !!sym(metric_col) < metric) %>%
        dplyr::pull(!!sym(name_col))
    out <- data.frame(features = c(up, down), direction = c(rep("up", length(up)), rep("down", length(down))))
    return(out)
}

#' Function to add missing rows to a matrix
#' Arguments:
#'   - df: matrix, rows = genes, cols = samples
#'   - rows: vector of rownames to add
#'   - sorted: logical, whether to sort the rows alphabetically
#' Returns:
#'   - df: matrix, rows = genes, cols = samples
add_missing_rows <- function(
    df, # cols = samples, rows = genes
    rows, # cols = samples, rows = genes
    sorted = TRUE ) {
  missingRowNames <-  rows[which(!rows %in% rownames(df))]
  print(missingRowNames)
  df_tmp <- as.data.frame(matrix(0,
                                 nrow = length(missingRowNames),
                                 ncol = dim(df)[2]
  )
  )
  # print(dim(df_tmp))
  # print(length(missingRowNames))
  colnames(df_tmp) <- colnames(df)
  rownames(df_tmp) <- missingRowNames
  # print(head(df_tmp))
  # print(head(df))
  df <- rbind(df_tmp, df)
  if (sorted) {
    df <- df[order(rownames(df)),]
  }
  
  return(df)
}

#' Function to make a SummarizedExperiment object
#' Arguments:
#'   - countsMatr: matrix of counts, rows = genes, cols = samples
#'   - colData: data.frame of sample metadata, rows = samples
#' Returns:
#'   - se: SummarizedExperiment object
make_se <- function(countsMatr, colData) {
  require(SummarizedExperiment)
  require(BiocGenerics)
  sample_ids <- intersect(colnames(countsMatr), row.names(colData))
  se <- SummarizedExperiment(assays = list(counts = as.matrix(countsMatr[,sample_ids])), 
                             colData = colData[sample_ids,])
  return(se)
}

#' Function to save a SummarizedExperiment object into multiple files if the slots are filled in the object
#' Arguments:
#'   - se: SummarizedExperiment object
#'   - path: path to save files
#'   - normalize: how to normalize the counts
save_se <- function(se, path, normalize = 'mor') {
  # make sure the path exists
  dir.create(path, showWarnings = F, recursive = T)

  # extract counts, colData, and rowData
  counts <- normalize_counts(se, method = normalize)
  colData <- as.data.frame(colData(se))
  rowData <- as.data.frame(rowData(se))

  # Save counts
  write.csv(counts, file = file.path(path, paste0('counts_', normalize,'.csv')), quote = F, row.names = T)
  # save the colData if it exists
  if (!is.null(colData)) {
    write.csv(colData, file = file.path(path, 'metadata.csv'), quote = T, row.names = T)
  }
  # save the rowData if it exists
  if (!is.null(rowData)) {
    write.csv(rowData, file = file.path(path, 'rowData.csv'), quote = T, row.names = T)
  }

}

#' Function to summarize a dataframe
#' Arguments:
#'   - df: dataframe to summarize
#' Returns:
#'   - df: dataframe with summary statistics
summarize_df <- function(df) {
  df %>%
    summary() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename('variable' = 'rowname') %>%
    mutate(variable = as.character(variable))

  return(df)
}

#' Function to normalize dds object and return counts
#' Arguments:
#'   - dds: DESeq2 object
#'   - method: normalization method
#' Returns:
#'   - counts: normalized counts
normalize_counts <- function(dds, method = 'mor', log2 = FALSE) {
  suppressPackageStartupMessages(require(DESeq2))
  suppressPackageStartupMessages(require(SummarizedExperiment))
  suppressPackageStartupMessages(require(BiocGenerics))
  suppressPackageStartupMessages(require(edgeR))
  suppressPackageStartupMessages(require(singscore))
  # Error handling
  options <- c('mor', 'vst', 'vsd', 'log2', 'rld', 'cpm', 'rlog', 'rpkm', 'none', 'tmm', 'log2-mor', 'rank')
  if (method %in% options) {
    normalize <- method
  } else {
    stop('Invalid normalization method. Please choose from: log2-mor, mor, vst, vsd, log2, rld, cpm, rlog, rpkm, tmm, rank, none')
  }

  # Get normalized counts
  if (normalize == "mor") {counttable <- counts(dds, normalize = T)}
  if (normalize %in% c("vsd", "vst")) {counttable <- assay(varianceStabilizingTransformation(dds))}
  if (normalize == "log2") {counttable <- log2(assay(dds)+1)}
  if (normalize == "log2-mor") {counttable <- log2(counts(dds, normalize = T)+1)}
  if (normalize == "rld") {counttable <- rlog(dds)}
  if (normalize == "cpm") {counttable <- cpm(dds)}
  if (normalize == "rlog") {counttable <- rlog(dds)}
  if (normalize == "rpkm") {counttable <- rpkm(dds)}
  if (normalize == "none") {counttable <- counts(dds)}
  if (normalize == "tmm") {counttable <- cpm(calcNormFactors(dds, method = "TMM"))}
  if (normalize == "rank") {counttable <- rankGenes(dds)}
  if (normalize == "none") {counttable <- counts(dds)}
  counts <- as.data.frame(counttable)
  if (log2) {
    counts <- log2(counts + 1)
  }
  return(counts)
}

#' Function to turn list of lists into a dataframe
#' Arguments:
#'   - list: list of lists
#' Returns:
#'   - df: dataframe
list_of_lists_to_df <- function(list) {
  df <- do.call(rbind, lapply(list, function(x) data.frame(x)))
  return(df)
}

# Function to remove NA variables from a summarized experiment object
# Arguments:
#   se: SummarizedExperiment object
#   columns: list of columns to remove NAs from
# Outputs:
#   DESeqDataSet object with NAs removed
remove_na_variables <- function(se, columns) {
  # Get the colData from the summarized experiment object
  col_data <- as.data.frame(colData(se))
  assay_data <- assay(se)
  
  # Update the DESeqDataSet object with the modified colData
  col_data <- drop_na(col_data, any_of(columns))
  col_data <- DataFrame(col_data)
  new_se <- make_se(assay_data, col_data)

  # Return the new DESeqDataSet object
  return(new_se)
}

#' Function to get the number of samples in each group
#' Arguments:
#'  - metadata: dataframe, metadata
#'  - id: character, column name of the sample ID
#'  - events_term: character, prefix of the events columns
#'  - subset: character, column name of the grouping variable
get_events_breakdown <- function(metadata, id = 'PATNUM', events_term = 'C_', subset = NULL) {
    # get the match of the PATNUMs
    metadata_patnums <- metadata[,id]
    if (!is.null(subset)) {
        patnum_match <- metadata_patnums[metadata[,subset] == 1]
    } else {
        patnum_match <- metadata_patnums
    }

    # catch if there are no matches
    if (is.na(patnum_match[1])) {
        stop('No matches found from the PATNUMs provided')
    }

    # Subset the metadata and break down the subtypes
    subtype_breakdown <- metadata[patnum_match,] %>% 
        select(starts_with(events_term)) %>%
        # count the 1s and 0s
        summarise_all(~ sum(., na.rm = TRUE)) %>%
        mutate( total = rowSums(.))
    return(subtype_breakdown)
}

#' Function to get make pairwise combinations
#' Arguments:
#'   - vec: vector of values
#' Returns:
#'   - list of pairwise combinations
pairwise_combos <- function(vec) {
    unique_vals <- unique(vec)
    combos <- combn(unique_vals, 2)
    list_of_combos <- list()
    for (i in 1:ncol(combos)) {
        list_of_combos[[i]] <- combos[, i]
    }
    return(list_of_combos)
}

#' Function to get the variable name
#' Arguments:
#'  - var: variable name
#' Returns:
#'  - varName: character, variable name
varName <- function(var) {
  deparse(substitute(var))
}

#' Function to one hot encode dataframe column
#' Arguments:
#'   - df: dataframe
#'   - column: character, column name
#'   - binary: logical, whether to make the column binary
#' Returns:
#'   - df: dataframe with one hot encoded column
one_hot_encode_ovr <- function(df, column, binary = TRUE) {
  # Get the unique values of the column
  unique_vals <- unique(df[[column]])
  
  # For each unique value, create a new binary column
  for (val in unique_vals) {
    # Create a new column name
    new_col_name <- paste0(column, "_", val)
    message(glue("Creating column: {new_col_name}"))
    
    # Create the new column in the data frame
    if (binary) {
      df[[new_col_name]] <- ifelse(df[[column]] == val, 1, 0)
    } else {
      df[[new_col_name]] <- factor(ifelse(df[[column]] == val, val, 'rest'), levels = c('rest', val))
    }
  }
  
  # Return the modified data frame
  return(df)
}


#' Function to install all required packages from a script
#' Arguments:
#'   - script: character, path to script
#'   - 
#' Returns:
#'   - none
install_packages_from_script <- function(script) {
  # Read in the script
  script <- readLines(script)
  
  # Get the packages
  packages <- script[grepl('library', script)]
  packages <- gsub('library\\(', '', packages)
  packages <- gsub('\\)', '', packages)
  packages <- gsub('"', '', packages)
  packages <- gsub("'", '', packages)
  packages <- gsub(' ', '', packages)

  # if packages is empty, try to get the packages from the require statements
  if (length(packages) == 0) {
    packages <- script[grepl('require', script)]
    packages <- gsub('require\\(', '', packages)
    packages <- gsub('\\)', '', packages)
    packages <- gsub('"', '', packages)
    packages <- gsub("'", '', packages)
    packages <- gsub(' ', '', packages)
  }

  # if that is still empty, try to get the packages from the pkg load statements
  if (length(packages) == 0) {
    packages <- script[grepl('package', script)]
    packages <- gsub('package\\(', '', packages)
    packages <- gsub('\\)', '', packages)
    packages <- gsub('"', '', packages)
    packages <- gsub("'", '', packages)
    packages <- gsub(' ', '', packages)
  }

  # try to install the packages and
  # if it fails, try to install with biocmanager
  # if that fails, print the error
  # capture the output
  output <- capture.output({
    for (pkg in packages) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
      }, error = function(e) {
        tryCatch({
          BiocManager::install(pkg, dependencies = TRUE)
        }, error = function(e) {
          print(e)
        })
      })
    }
  })
}