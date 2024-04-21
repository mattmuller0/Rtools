###########################################################################
#
#                            rna_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: rna_functions

#======================== LIBRARIES ========================
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(singscore))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R")

#======================== CODE ========================
#' Function to compare one column to many columns
#' Arguments:
#'   - df: data frame
#'   - col: column to compare
#'   - cols: columns to compare to
#'   - outdir: output directory
#'   - ...: additional arguments to pass to stats functions
#' Returns:
#'   - list of stats results and plots
compare_one_to_many <- function(df, col, cols, outdir, ...) {
    # set up output directory
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    data_type <- class(df[[col]])
    stats_list <- list()
    plot_list <- list()
    for (x in cols) {
        # check if column exists
        if (!(col %in% colnames(df))) {
            message(glue::glue("Column {col} not found in data frame."))
            next
        }

        base_plot <- ggplot(df, aes(x = !!sym(x), y = !!sym(col))) + labs(title = glue::glue("{col} vs {x}"), x = x, y = col)
        if (data_type == "numeric") {
            # numeric data
            stats_results <- df %>% rstatix::cor_test(!!sym(col), !!sym(cols))
            plot_results <- base_plot + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman")
        } else if (data_type %in% c("character", "factor")) {
            # get the number of unique values
            n_unique <- df %>% dplyr::pull(col) %>% unique() %>% length()
            if (n_unique == 2) {
                # categorical data with < 2 unique values
                stats_results <- df %>% rstatix::anova_test(!!sym(col), !!sym(cols))
                plot_results <- base_plot + geom_boxplot() + stat_compare_means(method = "t.test")
            } else {
                # categorical data with > 2 unique values
                stats_results <- df %>% rstatix::anova_test(!!sym(col), !!sym(cols))
                plot_results <- base_plot + geom_boxplot() + stat_compare_means(method = "anova")
            }
        } else {
            stop(glue::glue("Data type {data_type} not supported."))
        }
        # save results
        ggsave(glue::glue("{outdir}/{col}_vs_{x}.png"), plot_results)
        stats_list[[x]] <- stats_results
        plot_list[[x]] <- plot_results
        }
    return(list(stats = stats_list, plots = plot_list))
}
