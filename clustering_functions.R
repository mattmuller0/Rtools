# HEADER --------------------------------------------------
# Author: Matthew Muller
# 
# Date: 2023-01-11
# 
# Script Name: clustering functions
# 
# Notes:
# Lots of fun functions I've made to work with for clustering


# LOAD LIBRARIES ------------------------------------------
library(tidyverse)
library(glue)
library(ggplot2)
library(NMF)

# CODE BLOCK ----------------------------------------------
# Add code here

# Function to prep a counts matrix for nmf
# Arguments:
#   - counts_matr: matrix, rows = genes, cols = samples
#   - outfile: file to save output to
#   - ...: other arguments to pass to nmf
# Returns:
#   - output: list of nmf and random nmf results
#   - saves output to file
prep_nmf_data <- function(counts_matr) {
  # We are doing to log transform the data and then scale it.
  # To get around negative values, we are going to make a positive and negative matrix
  # and then combine them

  # First, log transform the data
  log_counts_matr <- log2(counts_matr + 1)

  # Next, scale the data
  scaled_counts_matr <- t(scale(t(log_counts_matr)))

  # Now, we need to make a positive and negative matrix
  # First, we need to get the dimensions of the matrix
  nrow <- dim(scaled_counts_matr)[1]
  ncol <- dim(scaled_counts_matr)[2]

  # Now, we can make the positive and negative matrices
  pos_matr <- matrix(0, nrow = nrow, ncol = ncol)
  neg_matr <- matrix(0, nrow = nrow, ncol = ncol)

  # Now, we can fill them
  pos_matr[scaled_counts_matr > 0] <- scaled_counts_matr[scaled_counts_matr > 0]
  neg_matr[scaled_counts_matr < 0] <- abs(scaled_counts_matr[scaled_counts_matr < 0])

  # Now, we can combine them
  combined_matr <- rbind(pos_matr, neg_matr)

  # Now, we can return the output
  return(combined_matr)
}


# Function to estimate the rank of a matrix
# Arguments:
#   - counts_matr: matrix, rows = genes, cols = samples
#   - outfile: file to save output to
#   - ranks: vector of ranks to test
#   - runs: number of runs to test
#   - options: options for nmf
#   - ...: other arguments to pass to nmf
# Returns:
#   - output: list of nmf and random nmf results
#   - saves output to file
nmf_estimator <- function(
    counts_matr,
    outfile,
    ranks = 2:8, runs = 30,
    options = 'v',
    ...
    ) {
  require(NMF)
  # Call nmf with some settings to make test the ranks chosen
  # We also are going to randomize the data to see how it compares to
  # data that is purely random in terms of residuals and other metrics
  # Ideally, the next step is to plot the NMF
  nmf_ranks_out <- nmf(counts_matr, ranks, nrun = runs, 
                       .opt = options, ...) 
  rng_ranks_out <- nmf(randomize(counts_matr), ranks, nrun = runs, 
                       .opt = options, ...) 

  output <- list('nmf_ranks_out' = nmf_ranks_out, 'rng_ranks_out' = rng_ranks_out)

  saveRDS(output, file = outfile)
  return(output)
}

# Function to plot the results of nmf_estimator
# Arguments:
#   - nmf_out: output from nmf_estimator
#   - outfile: file to save output to
#   - plot: whether to plot the results
# Returns:
#   - output: list of plots
#   - saves output to file
nmf_plotter <- function(nmf_out) {
  require(NMF)
  # Plot the results of nmf_estimator
  # We want to create a new metric that is the cophenetic correlation * dispersion
  # This will give us a better idea of how well the rank is doing

  # First, we need to get the cophenetic correlation
  cophenetic_correlation <- nmf_out$nmf_ranks_out$cophenetic
  # Next, we need to get the dispersion
  dispersion <- nmf_out$nmf_ranks_out$dispersion
  # Now, we can multiply them together
  cophenetic_dispersion <- cophenetic_correlation * dispersion

  # Now, we need to get the random cophenetic dispersion
  rng_cophenetic_correlation <- nmf_out$rng_ranks_out$cophenetic
  rng_dispersion <- nmf_out$rng_ranks_out$dispersion
  rng_cophenetic_dispersion <- rng_cophenetic_correlation * rng_dispersion

  # now combine them into a data frame
  cophenetic_df <- data.frame('rank' = nmf_out$nmf_ranks_out$rank,
                              'cophenetic' = cophenetic_correlation,
                              'dispersion' = dispersion,
                              'cophenetic_dispersion' = cophenetic_dispersion,
                              'type' = 'nmf')

  rng_cophenetic_df <- data.frame('rank' = nmf_out$rng_ranks_out$rank,
                                  'cophenetic' = rng_cophenetic_correlation,
                                  'dispersion' = rng_dispersion,
                                  'cophenetic_dispersion' = rng_cophenetic_dispersion,
                                  'type' = 'random')

  # combine them
  combined_df <- rbind(cophenetic_df, rng_cophenetic_df)

  # plot them
  cophenetic_plot <- ggplot(combined_df, aes(x = rank, y = cophenetic_dispersion, color = type)) +
    geom_point() + 
    labs(x = "Rank", y = "Cophenetic Dispersion") + 
    theme_matt()
}

# Function to estimate the kmeans clustering of a matrix using the elbow method
# Arguments:
#   - counts_matr: matrix, rows = genes, cols = samples
#   - outfile: file to save output to
#   - ks: vector of ks to test
#   - runs: number of runs to test
#   - ...: other arguments to pass to kmeans
# Returns:
#   - output: list of kmeans
#   - saves output to file
kmeans_estimator <- function(
    counts_matr,
    ks = 2:8,
    outfile,
    ...
    ) {
    km_list <- list()
    for (k in ks){
        km_list[[k]] <- kmeans(counts_matr, k, ...)
    }
    output <- list('km_list' = km_list)
    saveRDS(output, file = outfile)
    return(output)
}

# Function to estimate the hierarchical clustering of a matrix
# Arguments:
#   - counts_matr: matrix, rows = genes, cols = samples
#   - clusters: vector of clusters to test
#   - outfile: file to save output to
#   - ...: other arguments to pass to hclust
# Returns:
#   - output: list of clustering
#   - saves output to file
hclust_estimator <- function(
    counts_matr,
    clusters = 2:8,
    outfile,
    ...
    ) {
    hclust_list <- list()
    for (k in clusters){
        hclust_list[[k]] <- hclust(dist(counts_matr), ...)
    }
    output <- list('hclust_list' = hclust_list)
    saveRDS(output, file = outfile)
    return(output)
}

#======================== Old Code ========================#

# rank_estimator <- function(
#     expression_matrix,
#     ranks = c(2:8),
#     save_file = NULL
# ){
#   require(NMF)
#   # Runk NMF ranks
#   rank_out = nmf(expression_matrix, rank = ranks,
#                  method = nmf.getOption("default.algorithm"),
#                  nrun = 50)
  
#   # Graph cophenetic correlation plots
#   cophenetic_correlation_plot <- ggplot(rank_out, 
#                                         aes(x=rank, y=cophenetic)) + 
#     geom_point() + labs(x="Rank", y="Cophenetic Correlation Coefficient")
#   # Save data
#   if (!is.null(save_file)) {
#     save(rank_out, cophentic_correlation_plot,
#          file = save_file)
#   }
#   output <- list(rank_out, cophenetic_correlation_plot)
#   return(output)
# }
  
  
  
  
  
  



