###########################################################################
#
#                            rna_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: rna_functions

#======================== LIBRARIES ========================
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(singscore))
suppressMessages(library(tidyverse))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R")

message("
                                          .db--D,
                                          db( o\ o
                              _____       /  ,_   \
                       ,,---~~     ~~~~--;   /  ~-,$
                      /    __               /
                     /    /  \             /
                    (    ;    \            )
                     |  |       '--__;--   /
                     ) )                \ |
                    / /                 |:|
                   / /                  |;|
                   (||                  |:|
                   |,\                  |:!
 _________________.\\Le._______________.\\Le.__Pr59_____
")

#======================== Testing Functions ========================
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
    stats_list <- list()
    plot_list <- list()
    for (x in cols) {
        message(glue::glue("Comparing {col} to {x}..."))
        # check if column exists
        if (!(x %in% colnames(df))) {
            message(glue::glue("Column {col} not found in data frame."))
            next
        }

        data_type <- class(df[[x]])
        base_plot <- ggplot(tidyr::drop_na(df, dplyr::any_of(c(col, x))), aes(!!sym(x), !!sym(col))) + labs(title = glue::glue("{col} vs {x}"), x = x, y = col)
        if (data_type %in% c("integer", "numeric")) {
            # numeric data
            stats_results <- df %>%
                rstatix::cor_test(col, x, method = "spearman") %>%
                dplyr::select(variable = var2, method = method, estimate = cor, p = p) # nolint
            plot_results <- base_plot + geom_point() + geom_smooth(method = "lm") + stat_cor(method = "spearman")
        } else if (data_type %in% c("character", "factor")) {
            # get the number of unique values
            n_unique <- df %>% dplyr::pull(x) %>% na.omit() %>% unique() %>% length()
            if (n_unique == 2) {
                # categorical data with < 2 unique values
                stats_results <- df %>%
                    rstatix::t_test(as.formula(glue::glue("{col} ~ {x}"))) %>%
                    dplyr::mutate(variable = x, method = "t.test") %>%
                    dplyr::select(variable, method, estimate = statistic, p = p) # nolint
                plot_results <- base_plot + geom_boxplot() + stat_compare_means(method = "t.test")
            } else {
                # categorical data with > 2 unique values
                stats_results <- df %>%
                    rstatix::anova_test(as.formula(glue::glue("{col} ~ {x}"))) %>%
                    dplyr::as_tibble() %>%
                    dplyr::mutate(method = "anova") %>%
                    dplyr::select(variable = Effect, method, estimate = F, p = p) # nolint
                plot_results <- base_plot + geom_boxplot() + stat_compare_means(method = "anova")
            }
        } else {
            stop(glue::glue("Data type {data_type} not supported."))
        }
        # save results
        ggsave(glue::glue("{outdir}/{col}_vs_{x}.pdf"), plot_results)
        stats_list[[x]] <- stats_results
        plot_list[[x]] <- plot_results
        }
    stats_list <- dplyr::bind_rows(stats_list)
    write.csv(stats_list, glue::glue("{outdir}/stats_results.csv"))
    return(list(stats = stats_list, plots = plot_list))
}