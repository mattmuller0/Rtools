###########################################################################
#
#                            rna_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: rna_functions

#======================== LIBRARIES ========================#
suppressMessages(library(tidyverse))
suppressMessages(library(glue))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")

#======================== CODE ========================
# Function to take any codes and make a composite variable csv
#   Args:
#     lol: list of lists
#     composite_name: name of composite variable
#   Outputs:
#     csv of composite variable
make_compsite_coding_csv <- function(
    lol,
    composite_name
    ) {
    # ensure the lol args are all named vectors
    if (!all(sapply(list(lol), is.vector))) {
        stop("All unnamed arguments must be named vectors")
    }
    out <- data.frame(
        key = rep(composite_name, length(unlist(unlist(lol)))),
        field = c(rep(names(lol), sapply(lol, length))),
        value = c(unlist(unlist(lol)))
    )
    return(out)
}

# Function to make a phenotype from the composite encoding
#   Args:
#     hesin: data frame of hesin data pivoted longer
#            [dnx_hesin_id, Participant ID, field, value]
#     encoding: data frame of composite encoding [key, field, value]
#   Outputs:
#     data frame with composite encoding
make_phenotypes <- function(
    hesin,
    encoding
    ) {
    # ensure the required columns are present
    reqs <- c("Participant ID", "field", "value")
    if (!all(reqs %in% colnames(hesin))) {
        stop(glue::glue("Missing required columns in hesin data\n[{paste0(reqs, collapse = ", ")}]"))
    }
    reqs <- c("key", "field", "value")
    if (!all(reqs %in% colnames(encoding))) {
        stop(glue::glue("Missing required columns in encoding data\n[{paste0(reqs, collapse = ", ")}]"))
    }

    out <- list()
    for (i in seq_along(unique(encoding$key))) {
        k <- unique(encoding$key)[i]
        print(glue("Processing {unique(encoding$key)[i]}"))
        tmp <- encoding %>% filter(key == k)
        fields <- unique(tmp$field)
        o <- list()
        for (j in seq_along(fields)) {
            print(glue("    Field: {fields[j]}"))
            v <- unique(tmp %>% filter(field == fields[j]) %>% pull(value))
            v_reg <- paste0("^(", paste(v, collapse = "|"), ")")
            print(glue("        Values: {paste0(v, collapse = '|')}"))
            o[[j]] <- hesin %>%
                filter(field == fields[j]) %>%
                filter(grepl(v_reg, value)) %>%
                mutate(phenotype = k)
            print(glue("        N: {nrow(o[[j]])}"))
        }
        out[[i]] <- do.call(rbind, o)
    }
    out <- do.call(rbind, out)
    summary <- out %>%
        group_by(phenotype) %>%
        summarise(n = n())
    print("Summary of phenotypes")
    print(summary)
    return(out)
}