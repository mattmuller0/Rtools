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
library(tidyverse)

# LOAD FUNCTIONS
# space reserved for sourcing in functions

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

# generate a json file from a list
# Arguments:
#   - list: list, list to convert to json
#   - filename: character, name of the file to write
# Returns:
#   - none
list_to_json <- function(list, filename) {
  jsonlite::write_json(list, filename)
}