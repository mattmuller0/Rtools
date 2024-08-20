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
suppressMessages(library(tidyverse))

message("
            _.------.                        .----.__
           /         |_.       ._           /---.__  \
          |  O    O   |\\___  //|          /       `\ |
          |  .vvvvv.  | )   `(/ |         | o     o  |.
          /  |     |  |/      \ |  /|   ./| .vvvvv.   |\
         /   `^^^^^'  / _   _  `|_ ||  / /| |     |   | \
       ./  /|         | O)  O   ) \ | //' | `^vvvv'   |/\\
      /   / |         \        /  | | ~   \           |  \\
      \  /  |        / \ Y   /'   | \     |           |   ~
       `'   |  _     |  `._/' |   |  \     7          /
         _.-'-' `-'-'|  |`-._/   /    \ _ /     .    |
    __.-'            \  \   .   / |_.  \ -|_/\ / `--.|_
 --'                  \  \ |   /    |  |              `-
                       |uU |UU/     |  / 
                        `-.____.    |  | 
                                `--'  `-
")

# LOAD FUNCTIONS
# space reserved for sourcing in functions

###########################################################################
#
#                                 CODE
#
###########################################################################
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
  packages <- gsub('suppressMessages\\(', '', packages)
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
        install.packages(pkg, dependencies = TRUE, )
      }, error = function(e) {print(e)})
      tryCatch({
        BiocManager::install(pkg, dependencies = TRUE, update = TRUE, ask = FALSE)
      }, error = function(e) {print(e)})
    }
  })
  return(output)
}

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

# function to search a vector for a string
# Arguments:
#   - vector: vector, vector to search
#   - string: character, string to search for
# Returns:
#   - named vector of indices
search_vector <- function(vector, string, ...) {
  idx <- grep(string, vector, ...)
  vector <- vector[idx]
  names(idx) <- vector
  return(idx)
  return(vector)
}

drop_na_rows <- function(data, percent_allowed_missing) {
  if (percent_allowed_missing < 0 | percent_allowed_missing > 1) {
    stop("percent_allowed_missing must be between 0 and 1")
  }
  max_na_count <- percent_allowed_missing * ncol(data)
  data %>%
    rowwise() %>%
    filter(sum(is.na(c_across(everything()))) <= max_na_count) %>%
    ungroup()
}

drop_na_cols <- function(data, percent_allowed_missing) {
  if (percent_allowed_missing < 0 | percent_allowed_missing > 1) {
    stop("percent_allowed_missing must be between 0 and 1")
  }
  max_na_count <- percent_allowed_missing * nrow(data)
  data %>% select(where(~ sum(is.na(.)) <= max_na_count))
}