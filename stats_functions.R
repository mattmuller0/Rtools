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
# Add code here
#



# CORRELATIONS
# Function to calculate correlations between two SE objects
# Arguments:
#   - SE1: SummarizedExperiment object
#   - SE2: SummarizedExperiment object
#   - pval_filter: numeric, p-value threshold for filtering
#   - rval_filter: numeric, r-value threshold for filtering
#   - zval_filter: numeric, z-value threshold for filtering
#   - cor_method: character, correlation method to use
# Returns:
#   - correlations: data.frame, correlations between SE1 and SE2
sampleCorrelations <- function(
    SE1, SE2,
    pval_filter = 1,
    rval_filter = NA,
    zval_filter = NA,
    
    cor_method = 'pearson'
) {
  # SE1 <- SE1[rownames(SE1) %in% rownames(SE2), colnames(SE1) %in% colnames(SE2)]
  # SE2 <- SE2[rownames(SE2) %in% rownames(SE1), colnames(SE2) %in% colnames(SE1)]
  
  
  # This needs to be done in an apply loop of some sort for each row
  # also, make sure it's paired!!
  cor.test_ <- function(x) {cor.test(x, method=cor_method)}
  correlations <- mapply(cor.test, SE1, SE2)
  # correlations <- correlations[,correlations[3, ] <= pval_filter]
  # if (!is.na(rval_filter) & rval_filter >= 0) {correlations <- correlations[,correlations[4, ] >= rval_filter]}
  # if (!is.na(rval_filter) & rval_filter < 0) {correlations <- correlations[,correlations[4, ] < rval_filter]}
  # if (!is.na(zval_filter) & zval_filter >= 0) {correlations <- correlations[,correlations[1, ] >= zval_filter]}
  # if (!is.na(zval_filter) & zval_filter < 0) {correlations <- correlations[,correlations[1, ] < zval_filter]}
  # So this next bit is cause I am using mapply earlier, and clearly not processing data
  # as well as I could. This is a clear bandaid over the issue but hey it works lol
  rownames_ <- colnames(correlations)
  correlations <- as.data.frame(t(correlations))
  correlations <- sapply(correlations[-9], unlist) %>% as.data.frame()
  # Idk what is happening here but I'm doubling the number of rows
  # rownames(correlations) <- rownames_
  return(correlations)
}



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

  # Fix the table 1 by adding norm to the test column
  # table1 <- table1 %>% 
  #   mutate(
  #     test = case_when(
  #       p != '' & test != 'nonnorm' ~ 'norm',
  #       TRUE ~ test
  #     )
  #     )
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
    if (method %in% c('fisher', 'chisq')) {
        method <- method
    } else {
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

# Function to do a log rank test on a dataframe
# Arguments:
#   - data: data.frame, data to perform log rank test on
#   - comparisons: character, comparisons to perform log rank test on
#   - censors: character, censors to perform log rank test on
# Returns:
#   - p_values: data.frame, p-values from log rank test
log_rank_test <- function(
  data, 
  comparisons, 
  censors,

  censor_prefix = 'censor_',
  time_prefix = 'time_to_'
  ) {
  library(survival)
  
  # Create an empty dataframe to store the p-values
  p_values <- data.frame(matrix(ncol = length(comparisons), nrow = length(censors)))
  colnames(p_values) <- comparisons
  rownames(p_values) <- censors
  
  # Loop through each variable of interest
  for (censor in censors) {
    # Subset the data for the current variable of interest
    # get the corresponding time and event columns
    time <- gsub(censor_prefix, time_prefix, censor)

    # make sure censor and time are in the data
    if (!(time %in% colnames(data)) | !(censor %in% colnames(data))) {
      message(glue('time column {time} or censor column {censor} not in data'))
      next
    }
    
    # Loop through each comparison group
    for (i in 1:length(comparisons)) {
      # Subset the data for the current comparison group
      
      # Perform the log rank test
      # print(comparisons[i])
      fit <- survdiff(Surv(data[[time]], data[[censor]]) ~ data[[comparisons[i]]])
      p_value <- fit$p
      
      # Store the p-value in the p_values dataframe
      p_values[censor, comparisons[i]] <- p_value
    }
  }
  
  # Return the p-values dataframe
  return(p_values)
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
  # Subset the data
  # data <- data[, c(vars1, vars2)]
  
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

# Function to make survival curves and HR plots for a condition and censors
# Arguments:
#   - df: data.frame, data to make survival curves from
#   - condition: character, condition to make survival curves for
#   - censors: character, censors to make survival curves for
#   - outdir: character, output directory to save plots to
# Returns:
#   - survival_plots: list, list of survival plots
#   - HR_plots: list, list of HR plots
survival_analysis <- function(
  df, 
  condition, 
  censors, 
  outdir,

  image_type = 'png',
  censor_prefix = 'censor_',
  time_prefix = 'time_to_',

  survfit_args = list(),
  survdiff_args = list(),
  coxph_args = list()
  ) {
  require(ggsurvfit)
  require(survival)
  dir.create(file.path(outdir, 'survival_plots'), showWarnings = FALSE, recursive = TRUE)
  censors_basenames <- gsub(censor_prefix, '', censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  # verify the censors and time_vars are in the data
  # return the missing values in an error if not
  tryCatch({
    stopifnot(all(censors %in% colnames(df)))
    stopifnot(all(time_vars %in% colnames(df)))
  }, error = function(e) {
    message(glue('censor column {censors[!censors %in% colnames(df)]} not in data'))
    message(glue('time column {time_vars[!time_vars %in% colnames(df)]} not in data'))
    stop(e)
  })

  # make a list to store the survival plots
  survival_plots <- list()
  HR_list <- list()
  
  # loop over the censors
  for (idx in 1:length(censors)) {
    time_var <- time_vars[idx]
    censor <- censors[idx]
    cond_ <- df[[condition]]

    surv_obj <- survfit(Surv(df[[time_var]], df[[censor]]) ~ cond_)
    surv_diff <- survdiff(Surv(df[[time_var]], df[[censor]]) ~ cond_)
    coxmodel <- coxph(Surv(df[[time_var]], df[[censor]]) ~ cond_)

    surv_plot <- ggsurvfit(surv_obj) +
      coord_cartesian(clip = "off") +
      add_confidence_interval(alpha=0.1) +
      theme_classic(18) + 
      theme(legend.position = "top") +
      # add padding to axis for the risk table
      lims(x = c(0, max(surv_obj$time)+100)) +
      labs(
        x = 'Time (days)',
        title = glue('Survival Analysis for {censor} (n={sum(surv_obj$n)})')
        ) +
      # add text for coxmodel
      annotate(
        "text", 
        x = max(surv_obj$time)*0.75, 
        y = 1, 
        label = paste0('HR = ', signif(exp(coef(coxmodel)), 3), '; p = ', signif(summary(coxmodel)$coefficients[5], 1)), size = 7
      ) +
      # make the y axis title readable
      theme(axis.title.y = element_text(angle = 90, vjust = 0.5))

      # add the surv_plot to the list
      survival_plots[[censor]] <- surv_plot
      ggsave(filename = file.path(outdir, 'survival_plots', glue('survival_plot_{censor}.{image_type}')), plot = surv_plot)

      # make a dataframe of the OR & pvalue and HR & pvalues
      HR_list[[censor]] <- data.frame(
        logrank_chisq = signif(surv_diff$chisq, 4),
        logrank_pvalue = signif(surv_diff$pvalue, 4),
        HR = signif(exp(coef(coxmodel)), 4),
        HR_ci_lower = signif(exp(confint(coxmodel))[1], 4),
        HR_ci_upper = signif(exp(confint(coxmodel))[2], 4),
        HR_pvalue = signif(summary(coxmodel)$coefficients[5], 4)
      )
  }
  HR_df <- do.call(rbind, HR_list)
  write.csv(HR_df, file.path(outdir, 'HR_df.csv'))
  out <- list(
    survival_plots = survival_plots,
    HR_df = HR_df
  )
  return(out)
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

# Function to make HR tables for a condition and censors
# Arguments:
#   - df: data.frame, data to make survival curves from
#   - condition: character, condition to make survival curves for
#   - censors: character, censors to make survival curves for
#   - per_sd: logical, whether to standardize the condition
#   - ovr: logical, whether to do one vs rest
# Returns:
#   - HR_plots: list, list of HR plots
hazard_ratios_table <- function(
  df, 
  condition, 
  censors, 
  controls = NULL,

  per_sd = FALSE,
  ovr = FALSE,

  censor_prefix = 'C_',
  time_prefix = 'T_',

  survfit_args = list(),
  survdiff_args = list(),
  coxph_args = list()
  ) {
  require(ggsurvfit)
  require(survival)
  censors_basenames <- gsub(censor_prefix, '', censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  # verify the censors and time_vars are in the data
  # return the missing values in an error if not
  tryCatch({
    stopifnot(all(censors %in% colnames(df)))
    stopifnot(all(time_vars %in% colnames(df)))
  }, error = function(e) {
    message(glue('censor columns {paste0(censors[!censors %in% colnames(df)], collapse = ", ")} not in data\n'))
    message(glue('time column {paste0(time_vars[!time_vars %in% colnames(df)], collapse = ", ")} not in data\n'))
    stop(e)
  })

  # check if there are NA values in the condition
  if (any(is.na(df[, condition]))) {
    message('NA values in condition, removing')
    df <- df[!is.na(df[, condition]), ]
  }

  # Check if there are NA values in the controls
  if (!is.null(controls)) {
    if (any(is.na(df[, controls]))) {
      message('NA values in the controls')
    }
  }
  
  # Check if there are NA values in the censors
  if (any(is.na(df[, censors]))) {
    message('NA values in the censors')
  }

  # check if we are per_sd
  if (per_sd) {
    # make sure the condition is continuous
    if (!is.numeric(df[, condition])) {
      stop('condition must be numeric if per_sd is TRUE')
    }
    df[, condition ] <- scale(df[, condition ])
  }


  # make a list to store the survival plots
  HR_list <- list()
  
  if (!ovr) {
  # loop over the censors
    for (idx in 1:length(censors)) {
      time_var <- time_vars[idx]
      censor <- censors[idx]
      cond <- df[[condition]]
      if (!is.null(controls)) {
        # extract the controls from the df (there may be more than one)
        fmla <- as.formula(paste0('Surv(df[[time_var]], df[[censor]]) ~ ', condition, ' + ', paste(controls, collapse = ' + ')))
      } else {
        fmla <- as.formula(paste0('Surv(df[[time_var]], df[[censor]]) ~ ', condition))
      }

      surv_obj <- do.call(survfit, list(fmla, data = df))
      surv_diff <- do.call(survdiff, list(fmla, data = df))
      coxmodel <- do.call(coxph, list(fmla, data = df))

      tidy_cox <- broom::tidy(coxmodel)
      tidy_confint <- confint(coxmodel)
      tidy_surv <- broom::tidy(surv_obj)

      # make a dataframe of the OR & pvalue and HR & pvalues
      HR_list[[censor]] <- data.frame(
            censor = censor,
            condition = condition,
            HR = signif(exp(tidy_cox[1, 'estimate']), 4),
            HR_ci_lower = signif(exp(tidy_confint[1, 1]), 4),
            HR_ci_upper = signif(exp(tidy_confint[1, 2]), 4),
            HR_pvalue = signif(tidy_cox[1, 'p.value'], 4)
          )
    }
  } else {
    # one hot encode the condition
    message('one hot encoding condition')
    df_ <- one_hot_encode_ovr(df, condition, binary = FALSE)

    # get the names of the new columns
    vals <- unique(df_[, condition])
    columns_ovr <- colnames(df_)[colnames(df_) %in% paste0(condition, '_', vals)]

    # lapply the function over the columns_ovr and censors
    message('doing one vs rest')
    HR_list <- lapply(
      columns_ovr,
      function(x) {
        final <- data.frame()
        # loop over the censors
        for (idx in 1:length(censors)) {
          time_var <- time_vars[idx]
          censor <- censors[idx]
          cond <- df_[, x]

          if (!is.null(controls)) {
            # extract the controls from the df (there may be more than one)
            fmla <- as.formula(paste0('Surv(df[[time_var]], df[[censor]]) ~ ', x, ' + ', paste(controls, collapse = ' + ')))
          } else {
            fmla <- as.formula(paste0('Surv(df[[time_var]], df[[censor]]) ~ ', x))
          }

          surv_obj <- do.call(survfit, list(fmla, data = df_))
          surv_diff <- do.call(survdiff, list(fmla, data = df_))
          coxmodel <- do.call(coxph, list(fmla, data = df_))

          tidy_cox <- broom::tidy(coxmodel)
          tidy_confint <- confint(coxmodel)
          tidy_surv <- broom::tidy(surv_obj)

          # make a dataframe of the OR & pvalue and HR & pvalues
          out <- data.frame(
            censor = censor,
            condition = x,
            HR = unlist(signif(exp(tidy_cox[1, 'estimate']), 4)),
            HR_ci_lower = signif(exp(tidy_confint[1, 1]), 4),
            HR_ci_upper = signif(exp(tidy_confint[1, 2]), 4),
            HR_pvalue = signif(tidy_cox[1, 'p.value'], 4)
          )

          final <- rbind(final, out)
        }

        # return the HR_list
        return(final)
      }
    )
  }

  HR_df <- do.call(rbind, HR_list)

  # do some last minute cleaning
  rownames(HR_df) <- gsub('censor_', '', rownames(HR_df))
  rownames(HR_df) <- gsub('\\.', '__', rownames(HR_df))
  
  return(HR_df)
  }

# Function to a glm trend table
# Arguments:
#   - data: data.frame, data to make trend table from
#   - x: character, x variable
#   - ys: character, y variables
#   - verbose: logical, whether to print model summaries
#   - ...: other arguments to pass to glm
# Returns:
#   - out_dat: data.frame, trend table
trend_table <- function(data, x, ys, verbose = FALSE, ...) {
    # make an empty dataframe
    out_dat <- data.frame()
    for (y in ys) {
            # check if there are any NAs in the data
            # and if so, remove them
            # if (any(is.na(data[, c(x, y)]))) {
            #     tmpDat <- data %>% drop_na(any_of(c(x, y)))
            #     if (verbose) {
            #         print(paste0('Removed NAs from ', y))
            #         print(dim(data))
            #     }
            # }

            # make the model
            model <- glm(
                as.formula(paste0(y, ' ~ ', x)), 
                data = data %>% drop_na(any_of(c(x, y))), 
                ...
                )
            if (verbose) {
                print(glue('Model summary for {y} ~ {x}'))
                print(summary(model))
            }
            
            # tidy the model
            tidyModel <- broom::tidy(model) %>%
                mutate(
                    variable = y,
                    n = nrow(data %>% drop_na(any_of(c(x, y))))
                ) %>%
                filter(term == x) %>%
                dplyr::select(variable, n, estimate, std.error, statistic, p.value)
            # add to the out_dat
            out_dat <- rbind(out_dat, tidyModel)
    }
    return(out_dat)
}


# Function to make filtered hazard ratio tables
filtered_hazard_ratio_table <- function(
  data, 

  condition,
  risks, 
  censors, 
  
  censor_prefix = 'C_',
  time_prefix = 'T_',
  per_sd = TRUE,
  ovr = FALSE,

  verbose = FALSE, 
  ...
  ) {
  # make sure the risks are all characters
  if (!all(sapply(risks, is.character))) {
    stop('risks must be characters')
  }

  surv_risk_res <- map(
      risks, 
      function(x) {
          vals <- na.omit(unique(data[,x]))
          print(glue("Filtering for {x}"))
          res <- map(
              vals, 
              function(y) {
                  tmp <- data %>% filter(!!sym(x) == y)
                  if (verbose) {
                      message(glue("    Subfiltering {y}"))
                      message(glue("    N = {nrow(tmp)}"))
                  }
                  tryCatch({
                      out <- hazard_ratios_table(
                          df = tmp, 
                          condition = condition,
                          censors = censors,
                          censor_prefix = censor_prefix,
                          time_prefix = time_prefix,
                          per_sd = per_sd,
                          ovr = ovr,
                          ...
                          )
                      out$x <- x
                      out$y <- y
                      out$n <- nrow(tmp)
                      return(out)
                  },
                      error = function(e) {
                          if (verbose) {message(glue("    ERROR: {e}"))}
                          return(NULL)
                      }
                  )
              }
          )
          res <- do.call(rbind, res)
      }
  )
  surv_risk_res <- do.call(rbind, surv_risk_res)
  return(surv_risk_res)
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