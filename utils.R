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
# This is just some utility functions that I use in my R scripts. I put them
# here so I can source them in my scripts and not have to copy and paste them
# into each script.

###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "devtools",
  "BiocManager",
  "SummarizedExperiment"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T) # nolint
}

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')

###########################################################################
#
#                                 CODE
#
###########################################################################
# Function to set up a directory for output
# Arguments:
#   - outpath: path to output directory
#   - overwrite: logical, whether to overwrite the directory if it already exists
# Returns:
#   - outpath: path to output directory
setup_output_dir <- function(outpath, overwrite = F) {
  if (dir.exists(outpath)) {
    if (overwrite) {
      dir.create(outpath, showWarnings = F, recursive = T)
    } else {
      stop("Output directory already exists. Set overwrite = T to overwrite.")
    }
  } else {
    dir.create(outpath, showWarnings = F, recursive = T)
  }
  return(outpath)
}

# Function to set the seed for reproducibility
# Arguments:
#   - seed: integer, seed to set
# Returns:
#   - seed: integer, seed that was set
set_seed <- function(seed) {
  set.seed(seed)
  return(seed)
}

# Function to time a function
# Arguments:
#   - func: function to time
# Returns:
#   - time: time it took to run the function
time_function <- function(func) {
  start_time <- Sys.time()
  func
  end_time <- Sys.time()
  time <- end_time - start_time
  return(time)
}

# Function to detect the number of cores on a machine
# Arguments:
#   - none
# Returns:
#   - cores: integer, number of cores on the machine
detect_cores <- function() {
  cores <- parallel::detectCores()
  return(cores)
}

# Function to detect the number of threads on a machine
# Arguments:
#   - none
# Returns:
#   - threads: integer, number of threads on the machine
detect_threads <- function() {
  threads <- parallel::detectCores(logical = F)
  return(threads)
}