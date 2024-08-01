# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2023-01-15
# 
# Script Name: Converting Functions Script
# 

# LOAD LIBRARIES ------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
suppressMessages(library(ggplot2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
# suppressMessages(library(AnnotationDbi))

#======================== ID Conversions ========================#
# Function to map gene IDs with annotationDbi
# Arguments:
#   - geneList: list or vector of gene IDs
#   - from: character, type of gene ID to convert from (if NULL, will detect)
#   - to: character, type of gene ID to convert to
#   - keepUnmatched: logical, keep unmatched gene IDs
#   - ...: additional arguments to pass to mapIds
# Returns:
#   - out_list: list, converted gene IDs
map_gene_ids  <- function(geneList, from = NULL, to, orgDb = org.Hs.eg.db, remove_missing = FALSE, ...) {
    # require("AnnotationDbi")
    require("org.Hs.eg.db")
    # if from is NULL, detect the type of gene ID
    if (is.null(from)) {
        from <- detect_gene_id_type(geneList)
    }

    # copy the gene list
    geneList_copy <- geneList

    # check if the IDs have . in them (implies versioning)
    if (grepl("\\.", geneList[1])) {
        geneList <- sapply(geneList, function(x) str_split(x, "\\.", simplify = TRUE)[1])
    }

    # convert gene ids
    out_list <- AnnotationDbi::mapIds(orgDb, keys = geneList, column = to, keytype = from, ...)

    # ake any NA values into the original gene id
    if (remove_missing) {
        out_list <- out_list[!is.na(out_list)]
    } else {
        out_list[is.na(out_list)] <- geneList_copy[is.na(out_list)] %>% as.vector()
    }

    # return the converted gene ids
    return(out_list)
    }

# Function to convert gene IDs with biomart
# Arguments:
#   - geneList: list or vector of gene IDs
#   - from: character, type of gene ID to convert from (defaults to ensembl_gene_id)
#   - to: character, type of gene ID to convert to (defaults to external_gene_name)
#   - uniqueRows: logical, keep only unique rows
#   - mart: biomaRt object, biomart to use (defaults to ensembl)
#   - ...: additional arguments to pass to getBM
# Returns:
#   - data.frame, converted gene IDs
convert_gene_ids_df <- function(
    geneList, 
    from = 'ensembl_gene_id', to = 'external_gene_name', 
    uniqueRows = TRUE, 
    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"),
    ...) {
  require("biomaRt")

  # copy the gene list
  geneList_copy <- geneList

  # check if the IDs have . in them (implies versioning)
  if (grepl("\\.", geneList[1])) {
    geneList <- sapply(geneList, function(x) str_split(x, "\\.", simplify = TRUE)[1])
  }

  # convert gene ids
  out_frame <- getBM(
    attributes = c(from, to), 
    filters = from, 
    values = geneList, 
    mart = mart, 
    uniqueRows = uniqueRows, 
    ...
    )

  return(out_frame)
}

# Function to detect gene ID type based on first gene ID
# Arguments:
#   - geneList: list or vector of gene IDs
# Returns:
#   - type: character, type of gene ID
detect_gene_id_type <- function(geneList, strip = TRUE) {
  # ID types regex:
  id_types <- list(
    ENSEMBL = "^ENSG[0-9]+$",
    ENSEMBLTRANS = "^ENST[0-9]+$",
    ENSEMBLPROT = "^ENSP[0-9]+$",
    ENTREZID = "^[0-9]+$",
    IMAGE = "^IMAGE:[0-9]+$",
    GOID = "^GO:[0-9]+$",
    PFAM = "^PF[0-9]+$",
    REFSEQ = "^N[MP]_[0-9]+$",
    ENZYME = "^[0-9]+(\\.(([0-9]+)|-)+)3$",
    MAP = "^[0-9XY]+((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9XY]+)?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$",
    UNIGENE = "^Hs\\.([0-9]+)$",
    SYMBOL = "^[A-Z][A-Z0-9]*(_[A-Z0-9]+)*(-[A-Z0-9]+)*$",
    GENEBANK_Nucleotide = "^[A-Z][0-9]5$",
    GENEBANK_Protein = "^[A-Z]3[0-9]5$",
    GENEBANK_WGS = "^[A-Z]4[0-9]8[0-9]?[0-9]?$",
    GENEBANK_MGA = "^[A-Z]5[0-9]7$",
    GENENAME = " ",
    .Affymetrix = "(^AFFX-)|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$",
    .Illumina = "^ILMN_[0-9]+$",
    .Agilent = "^A_[0-9]+_P[0-9]+$"
    )

  # check if any gene ids has a . in it and strip if so
  if (grepl("\\.", geneList[1]) && strip) {
    geneList <- sapply(geneList, function(x) str_split(x, "\\.", simplify = TRUE)[1])
    message("Gene ID version numbers have been stripped.")
  }

  # loop through the id types and get the matches for each gene
  matches <- sapply(id_types, function(x) grepl(x, geneList))

  # get the colnames when there is a match
  type <- colnames(matches)[which(matches[1,] == TRUE)] # nolint

  # if there are multiple matches, give a warning
  if (length(type) > 1) {
    warning(glue("Multiple gene ID types detected. Returning first match: {type[1]}"))
    type <- type[1]
  }

  # stop if no matches
  if (length(type) == 0) {
    stop("No gene ID type detected.")
  }

  # return the type
  return(type)
}

# Function to convert gene IDs
# Arguments:
#   - df: data.frame, data frame to subset
#   - type: character, type of gene ID to convert to
#   - keys: character, type of gene ID to convert from
# Returns:
#   - df2: data.frame, subsetted data frame
getMatrixWithSelectedIds <- function(df, type, keys) {
  # require("AnnotationDbi")
  require("org.Hs.eg.db")
  geneSymbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=rownames(df), column=c(type), keytype=keys, multiVals="first")
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  # subset your data frame based on the found_genes
  df2 <- df[names(found_genes), ]
  rownames(df2) <- found_genes
  return(df2)
}

#======================== Normalizations ========================#
# Function for mor normalizing data
# Arguments:
#   - data: data.frame, data to normalize
# Returns:
#   - manually_normalized: data.frame, normalized data
mor_normalization <- function(data) {
  require(dplyr)
  require(tibble)
  
  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    dplyr::filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}

# Function to convert RPKM to counts
# Arguments:
#   - rpkm_data: data.frame, data to convert
# Returns:
#   - counts_data: data.frame, converted data
rpkm_to_counts <- function(rpkm_data) {
  # Get the scaling factor for each sample
  scaling_factors <- 1e6 / colSums(rpkm_data)
  
  # Convert RPKM to counts
  counts_data <- round(rpkm_data * scaling_factors)
  
  # Return the counts data
  return(counts_data)
}
