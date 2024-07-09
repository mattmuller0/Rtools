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

#======================== Data Preprocessing Functions ========================
#' Function to select top genes by variance
#' Arguments:
#'  - df: data frame [samples x genes]
#'  - n: number of top genes to select
#' Returns:
#'  - data frame with columns [genes, variance]
select_top_variability <- function(df, n = 1000) {
    variances <- apply(df, 2, var)
    top_genes <- names(sort(variances, decreasing = TRUE)[1:n])
    return(data.frame(genes = top_genes, variance = variances[top_genes]))
}

#' Function to filter genes by expression
#' Arguments:
#' - df: data frame [samples x genes]
#' - min_expr: minimum expression value
#' Returns:
#' - data frame with filtered genes [genes, expression]
select_top_expression <- function(df, min_expr = 1) {
    exprs <- apply(df, 2, function(x) sum(x > min_expr))
    filtered_genes <- names(exprs[exprs > 0])
    return(data.frame(genes = filtered_genes, expression = exprs[filtered_genes]))
}

#' Function to filter genes by lasso regression
#' Arguments:
#' - df: data frame [samples x genes]
#' - y: response variable
#' - lambda: lambda value for lasso regression
#' - nfolds: number of cross-validation folds
#' - ...: additional arguments to pass to lasso function
#' Returns:
#' - data frame with filtered genes [genes, coefficient]
#' - data frame with lasso results
select_top_lasso <- function(df, y, lambda = NULL, nfolds = 10, ...) {
    requireNamespace("glmnet", quietly = TRUE)
    lasso_res <- glmnet::cv.glmnet(df, y, alpha = alpha, lambda = lambda, nfolds = nfolds, ...)
    message(glue::glue("Selected lambda: {lasso_res$lambda.min}"))
    coef <- as.data.frame(coef(lasso_res, s = "lambda.min"))
    coef <- coef[coef[, 1] != 0, ]
    return(coef)
}



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

#======================== Eigengene Functions ========================
#' Function to calculate eigengenes by principal component analysis
#' Arguments:
#'  - df: data frame [samples x genes]
#'  - outdir: output directory
#'  - ...: additional arguments to pass to stats functions
#' Returns:
#' - dataframe with eigengenes
eigengenes_pca <- function(df, outdir, pcs = 1:3, align_avg_expr = FALSE, ...) {
    requireNamespace("ggbiplot", quietly = TRUE)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # run PCA
    pca_res <- prcomp(df, ...)

    biplot <- ggbiplot::ggbiplot(pca_res, obs.scale = 1, var.scale = 0.5, groups = NULL, ellipse = TRUE)
    ggsave(glue::glue("{outdir}/biplot.pdf"), biplot)
    
    eigengenes <- as.data.frame(pca_res$x[, pcs])
    colnames(eigengenes) <- glue::glue("PC{pcs}")

    # align average expression
    if (align_avg_expr) {
        avg_expr <- rowMeans(df)
        corrs <- cor(eigengenes, avg_expr)
        print(corrs)
        eigengenes <- as.vector(sign(corrs)) * eigengenes
    }

    write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Function to calculate eigengenes by singular value decomposition (WIP)
#' Arguments:
#' - df: data frame [samples x genes]
#' - outdir: output directory
#' - pcs: number of principal components to return
#' - ...: additional arguments to pass to stats functions
#' Returns:
#' - dataframe with eigengenes
eigengenes_svd <- function(df, outdir, pcs = 1:3, align_avg_expr = FALSE, ...) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # run SVD
    svd_res <- svd(df, ...)

    eigengenes <- as.data.frame(svd_res$v[, pcs])
    colnames(eigengenes) <- glue::glue("PC{pcs}")

    # align average expression
    if (align_avg_expr) {
        avg_expr <- rowMeans(df)
        corrs <- cor(eigengenes, avg_expr)
        print(corrs)
        eigengenes <- as.vector(sign(corrs)) * eigengenes
    }

    write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Function to calculate eigengenes by non-negative matrix factorization (WIP)
#' Arguments:
#' - df: data frame [samples x genes]
#' - outdir: output directory
#' - pcs: number of principal components to return
#' - ...: additional arguments to pass to stats functions
#' Returns:
#' - dataframe with eigengenes
eigengenes_nmf <- function(df, outdir, pcs = 1:3, align_avg_expr = FALSE, ...) {
    requireNamespace("NMF", quietly = TRUE)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # run NMF
    nmf_res <- NMF::nmf(as.matrix(df), rank = pcs, seed = 420)

    eigengenes <- as.data.frame(NMF::basis(nmf_res))
    colnames(eigengenes) <- glue::glue("PC{pcs}")

    # align average expression
    if (align_avg_expr) {
        avg_expr <- rowMeans(df)
        corrs <- cor(eigengenes, avg_expr)
        print(corrs)
        eigengenes <- as.vector(sign(corrs)) * eigengenes
    }

    write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Function to calculate eigengenes by independent component analysis (WIP)
#' Arguments:
#' - df: data frame [samples x genes]
#' - outdir: output directory
#' - pcs: number of principal components to return
#' - ...: additional arguments to pass to stats functions
#' Returns:
#' - dataframe with eigengenes
eigengenes_ica <- function(df, outdir, pcs = 1:3, align_avg_expr = FALSE, ...) {
    requireNamespace("fastICA", quietly = TRUE)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # run ICA
    ica_res <- fastICA::fastICA(df, n.comp = pcs, ...)

    eigengenes <- as.data.frame(ica_res$S)
    colnames(eigengenes) <- glue::glue("PC{pcs}")

    # align average expression
    if (align_avg_expr) {
        avg_expr <- rowMeans(df)
        corrs <- cor(eigengenes, avg_expr)
        print(corrs)
        eigengenes <- as.vector(sign(corrs)) * eigengenes
    }

    write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}

#' Function to calculate eigengenes by a generalized linear model
#' Arguments:
#' - df: data frame [samples x genes]
#' - response: response variable
#' - outdir: output directory
#' - scale: logical, scale the data
#' - center: logical, center the data
#' - ...: additional arguments to pass to stats functions
#' Returns:
#' - dataframe with eigengenes
eigengenes_glm <- function(df, response, outdir, family = "gaussian", type.measure = "mse", ...) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # run GLM
    glm_res <- glmnet::glmnet(df, response, family = family, type.measure = type.measure, ...)
    glm_tidy <- broom::tidy(glm_res)
    saveRDS(glm_res, glue::glue("{outdir}/glmnet_model.rds"))
    
    pdf(glue::glue("{outdir}/glmnet_plot.pdf"))
    plot(glm_res)
    dev.off()

    # get the eigengenes as the predictions from the model
    eigengenes <- predict(glm_res, newx = df)
    eigengenes <- as.data.frame(eigengenes)
    print(eigengenes)
    write.csv(eigengenes, glue::glue("{outdir}/eigengenes.csv"))
    return(eigengenes)
}