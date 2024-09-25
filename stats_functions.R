###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-03-10
# 
# Script Name: Statistics Functions
# 
# Notes:
# general stats functions I use for stuff. most use SEs or DFs

###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(tibble))
suppressMessages(library(tableone))
suppressMessages(library(glue))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')

###########################################################################
#
#                                 CODE
#
###########################################################################
# Function to make a table 1
# Arguments:
#   - data: data.frame, data to make table 1 from
#   - groups: character, groups to stratify by
#   - vars: character, variables to include in table 1
#   - printArgs: list, arguments to pass to print
#   - ...: other arguments to pass to CreateTableOne
# Returns:
#   - table1: data.frame, table 1
stats_table <- function(
    data,
    groups,
    vars = NULL,
    printArgs = list(
      showAllLevels = FALSE,
      printToggle = FALSE
    ),
    ...
) {
  require(tableone)
  
  # Get all the the variables if none are specified
  if (is.null(vars)) {
    vars <- colnames(data)
    vars <- vars[!vars %in% groups]
  }

  # Make a table 1
  table1 <- CreateTableOne(
    vars = vars,
    strata = groups,
    data = data,
    addOverall = TRUE,
    ...
  )
  
  # Make a nice table 1
  table1 <- do.call(print, c(list(table1), printArgs))
  return(table1)
}

# Function to do an odds ratio between a gene set and hypergeometric enrichment results
# Arguments:
#   - gse: gse object from hypergeometric enrichment analysis in clusterProfiler
#   - ...: other arguments to pass to oddsratio
# Returns:
#   - odds_ratio: data.frame, enrichment dataframe with odds ratio and p-value added
hypergeometric_scoring <- function(
    gse,
    method = 'fisher',
    ...
    ) {
    # Get the geneset
    geneset <- gse@gene

    # Get the enrichment sets
    enrichment <- gse@geneSets

    # get the universe genes
    universe <- gse@universe

    # make a dataframe where each row is the universe
    df <- data.frame(
        row.names = universe,
        in_geneset = universe %in% geneset
        )
    # now add each enrichment set as a new column
        for (i in 1:length(enrichment)) {
            df[, names(enrichment)[i]] <- universe %in% enrichment[[i]]
        }

    # error handling on the method
    if (!method %in% c('fisher', 'chisq')) {
        stop('method must be either fisher or chisq')
    }
    
    # let's do the OR test for each enrichment set
    if (method == 'fisher') {
        ORs <- sapply(
            df[, -1], 
            function(x) fisher.test(table(x, df[, 'in_geneset']), ...)$estimate
            )
    } else if (method == 'chisq') {
        ORs <- sapply(
            df[, -1], 
            function(x) chisq.test(table(x, df[, 'in_geneset']), ...)$estimate
            )
    }

    # fix the names
    names(ORs) <- gsub('\\..*', '', names(ORs))

    # add the ORs to the enrichment results
    results_pathways <- rownames(gse@result)
    gse@result$odds_ratio <- ORs[results_pathways]

    # return the gse object
    return(gse)
}

# Function to make a correlation matrix of given data and variables
# Arguments:
#   - data: data.frame, data to make correlation matrix from
#   - vars1: character, variables to correlate
#   - vars2: character, variables to correlate
#   - method: character, correlation method to use
# Returns:
#   - cor_mat: data.frame, correlation matrix
correlation_matrix <- function(data, vars1, vars2, method = 'pearson', use = 'pairwise.complete.obs', ...) {
  # make a placeholder for the correlation matrix
  cor_mat <- data.frame()
  p_mat <- data.frame()
  # map the correlation function over the data
  cor_res <- map(
    vars1,
    function(x) {
      map(
        vars2,
        function(y) {
          res <- cor.test(data[, x], data[, y], method = method, use = use, ...)
          cor_mat[x, y] <<- res$estimate
          p_mat[x, y] <<- res$p.value
        }
      )
    }
  )
  # make a correlation matrix filtered by p-value
  cor_mat_filtr <- cor_mat
  cor_mat_filtr[p_mat > 0.05] <- NA
  return(list(cor_matr = cor_mat, pvalue_matr = p_mat, cor_matr_filtr = cor_mat_filtr))
}

# Function to softmax a vector
# Arguments:
#   - x: numeric, vector to softmax
# Returns:
#   - softmax: numeric, softmaxed vector
softmax <- function(x) {
    exp(x) / sum(exp(x))
}

# Function to min-max normalize a vector
# Arguments:
#   - x: numeric, vector to min-max normalize
# Returns:
#   - min_max_norm: numeric, min-max normalized vector
min_max_norm <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

#======================== Old Code ========================#
# # Function to calculate p-values from a correlation matrix
# # Arguments:
# #   - cor_matrix: data.frame, correlation matrix
# #   - n: numeric, number of samples
# # Returns:
# #   - p_values: data.frame, p-values from correlation matrix
# correlation_pvalues <- function(cor_matrix, n) {
#   t_values <- cor_matrix * sqrt((n - 2)/(1 - cor_matrix^2))
#   p_values <- 2 * pt(abs(t_values), df = n - 2, lower.tail = FALSE)
#   p_values[is.na(p_values)] <- 1
#   return(p_values)
# }

# # Function to make a correlation matrix of given data and variables
# # Arguments:
# #   - data: data.frame, data to make correlation matrix from
# #   - vars1: character, variables to correlate
# #   - vars2: character, variables to correlate
# #   - method: character, correlation method to use
# # Returns:
# #   - cor_mat: data.frame, correlation matrix
# correlation_matrix <- function(data, vars1, vars2, method = 'pearson', use = 'pairwise.complete.obs') {
#   # Subset the data
#   data <- data[, c(vars1, vars2)]
  
#   # Calculate the correlation matrix
#   cor_mat <- cor(
#     data[, vars1],
#     data[, vars2],
#     method = method,
#     use = use
#     )

#   # Calculate the p-values
#   # p_mat <- correlation_pvalues(cor_mat, nrow(data))
  
#   # Return the correlation matrix
#   return(cor_mat)
# }


# Wilcoxan Ranked Sum testing on genes in two summarized experiments
# gene_wilcox_test <- function(dds, design) {
#   # Extract count data from the DESeq object
#   count_data <- assay(dds)
  
#   # Extract metadata from the DESeq object
#   meta <- colData(dds) %>% as.data.frame() %>% dplyr::select(sym(condition), age, race, sex)
  
#   # Create a design matrix
#   design_matrix <- model.matrix( ~ age + sex + race, data = meta)
  
#   # Perform Wilcoxon ranked sum test while adjusting for age, sex, and race
#   fit <- apply(count_data, 1, function(x) {
#     model <- glm(x ~ design_matrix)
#     wilcox <- wilcox.test(model$residuals ~ meta[, condition], exact = F)
#     return(wilcox)
#   })  
  
#   # Extract p-values and adjust for multiple testing using the Benjamini-Hochberg method
#   res <- data.frame(pvalue = sapply(fit, function(x) c(x$p.value))) 
#   res <- res %>% mutate(padj = p.adjust(res$pvalue, method="BH"))
  
#   # View results
#   return(res)
# }

# # Function to a glm trend table
# # Arguments:
# #   - data: data.frame, data to make trend table from
# #   - x: character, x variable
# #   - ys: character, y variables
# #   - verbose: logical, whether to print model summaries
# #   - ...: other arguments to pass to glm
# # Returns:
# #   - out_dat: data.frame, trend table
# trend_table <- function(data, x, ys, verbose = FALSE, ...) {
#     # make an empty dataframe
#     out_dat <- data.frame()
#     for (y in ys) {
#             # check if there are any NAs in the data
#             # and if so, remove them
#             # if (any(is.na(data[, c(x, y)]))) {
#             #     tmpDat <- data %>% drop_na(any_of(c(x, y)))
#             #     if (verbose) {
#             #         print(paste0('Removed NAs from ', y))
#             #         print(dim(data))
#             #     }
#             # }

#             # make the model
#             model <- glm(
#                 as.formula(paste0(y, ' ~ ', x)), 
#                 data = data %>% drop_na(any_of(c(x, y))), 
#                 ...
#                 )
#             if (verbose) {
#                 print(glue('Model summary for {y} ~ {x}'))
#                 print(summary(model))
#             }
            
#             # tidy the model
#             tidyModel <- broom::tidy(model) %>%
#                 mutate(
#                     variable = y,
#                     n = nrow(data %>% drop_na(any_of(c(x, y))))
#                 ) %>%
#                 filter(term == x) %>%
#                 dplyr::select(variable, n, estimate, std.error, statistic, p.value)
#             # add to the out_dat
#             out_dat <- rbind(out_dat, tidyModel)
#     }
#     return(out_dat)
# }

# # CORRELATIONS
# # Function to calculate correlations between two SE objects
# # Arguments:
# #   - SE1: SummarizedExperiment object
# #   - SE2: SummarizedExperiment object
# #   - pval_filter: numeric, p-value threshold for filtering
# #   - rval_filter: numeric, r-value threshold for filtering
# #   - zval_filter: numeric, z-value threshold for filtering
# #   - cor_method: character, correlation method to use
# # Returns:
# #   - correlations: data.frame, correlations between SE1 and SE2
# sampleCorrelations <- function(
#     SE1, SE2,
#     pval_filter = 1,
#     rval_filter = NA,
#     zval_filter = NA,
    
#     cor_method = 'pearson'
# ) {
#   # SE1 <- SE1[rownames(SE1) %in% rownames(SE2), colnames(SE1) %in% colnames(SE2)]
#   # SE2 <- SE2[rownames(SE2) %in% rownames(SE1), colnames(SE2) %in% colnames(SE1)]
  
#   # This needs to be done in an apply loop of some sort for each row
#   # also, make sure it's paired!!
#   cor.test_ <- function(x) {cor.test(x, method=cor_method)}
#   correlations <- mapply(cor.test, SE1, SE2)
#   # correlations <- correlations[,correlations[3, ] <= pval_filter]
#   # if (!is.na(rval_filter) & rval_filter >= 0) {correlations <- correlations[,correlations[4, ] >= rval_filter]}
#   # if (!is.na(rval_filter) & rval_filter < 0) {correlations <- correlations[,correlations[4, ] < rval_filter]}
#   # if (!is.na(zval_filter) & zval_filter >= 0) {correlations <- correlations[,correlations[1, ] >= zval_filter]}
#   # if (!is.na(zval_filter) & zval_filter < 0) {correlations <- correlations[,correlations[1, ] < zval_filter]}
#   # So this next bit is cause I am using mapply earlier, and clearly not processing data
#   # as well as I could. This is a clear bandaid over the issue but hey it works lol
#   rownames_ <- colnames(correlations)
#   correlations <- as.data.frame(t(correlations))
#   correlations <- sapply(correlations[-9], unlist) %>% as.data.frame()
#   # Idk what is happening here but I'm doubling the number of rows
#   # rownames(correlations) <- rownames_
#   return(correlations)
# }
