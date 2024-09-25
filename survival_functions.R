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

  censor_prefix = 'C_',
  time_prefix = 'T_'
  ) {
  require(survival)
  
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
      message(glue('Missing either time or censor column for {censor}'))
      next
    }
    
    # Loop through each comparison group
    for (i in seq_len(comparisons)) {
      # Subset the data for the current comparison group
      comparison <- comparisons[i]
      # Perform the log rank test
      fmla <- as.formula(glue("Surv({time}, {censor}) ~ {comparison}"))
      fit <- do.call(survdiff, fmla, data = data)
      p_value <- fit$p
      
      # Store the p-value in the p_values dataframe
      p_values[censor, comparison] <- p_value
    }
  }
  
  # Return the p-values dataframe
  return(p_values)
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

    image_type = 'pdf',
    censor_prefix = 'C_',
    time_prefix = 'T_'
    ) {
    require(ggsurvfit)
    require(survival)
    dir.create(file.path(outdir, 'survival_plots'), showWarnings = FALSE, recursive = TRUE)
    censors_basenames <- gsub(censor_prefix, '', censors)
    time_vars <- paste0(time_prefix, censors_basenames)

    # verify the censors and time_vars are in the data
    # return the missing values in an error if not
    stopifnot(all(censors %in% colnames(df)))
    stopifnot(all(time_vars %in% colnames(df)))

    # make a list to store the survival plots
    survival_plots <- list()
    HR_list <- list()

    # loop over the censors
    for (idx in 1:length(censors)) {
    time_var <- time_vars[idx]
    censor <- censors[idx]
    fmla <- as.formula(glue("Surv({time_var, {censor}} ~ condition)"))

    surv_obj <- do.call(survfit, fmla, data = df)
    surv_diff <- do.call(survdiff, fmla, data = df)
    coxmodel <- do.call(coxph, fmla, data = df)

    surv_plot <- ggsurvfit(surv_obj) +
        coord_cartesian(clip = "off") +
        add_confidence_interval(alpha=0.1) +
        theme_classic(18) + 
        theme(legend.position = "top") +
        lims(x = c(0, max(surv_obj$time)+100)) +
        labs(x = 'Time (days)', title = glue('Survival Analysis for {censor} (n={sum(surv_obj$n)})')) +
        annotate("text", x = max(surv_obj$time)*0.75, y = 1, label = paste0('HR = ', signif(exp(coef(coxmodel)), 3), '; p = ', signif(summary(coxmodel)$coefficients[5], 1)), size = 7) +
        theme(axis.title.y = element_text(angle = 90, vjust = 0.5))

        # add the surv_plot to the list
        survival_plots[[censor]] <- surv_plot
        ggsave(filename = file.path(outdir, 'survival_plots', glue('survival_plot_{censor}.{image_type}')), plot = surv_plot)

        # make a dataframe of the OR & pvalue and HR & pvalues
        o <- data.frame(
        logrank_chisq = signif(surv_diff$chisq, 4),
        logrank_pvalue = signif(surv_diff$pvalue, 4),
        hazard_ratio = signif(exp(coef(coxmodel)), 4),
        ci_lower = signif(exp(confint(coxmodel))[1], 4),
        ci_upper = signif(exp(confint(coxmodel))[2], 4),
        hr_pvalue = signif(summary(coxmodel)$coefficients[5], 4)
        )
        HR_list[[censor]] <- o
    }
    HR_df <- do.call(rbind, HR_list)
    write.csv(HR_df, file.path(outdir, 'hazard_ratios.csv'))
    out <- list(survival_plots = survival_plots, hazard_ratios = HR_df)
    return(out)
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
  ) {
  require(ggsurvfit)
  require(survival)
  censors_basenames <- gsub(censor_prefix, '', censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  # verify the censors and time_vars are in the data
  # return the missing values in an error if not
  stopifnot(all(censors %in% colnames(df)))
  stopifnot(all(time_vars %in% colnames(df)))

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
    df[, condition] <- scale(df[, condition])
  }

  # make a list to store the survival plots
  HR_list <- list()
  if (!ovr) {
  # loop over the censors
    for (idx in seq_len(censors)) {
      time_var <- time_vars[idx]
      censor <- censors[idx]
      cond <- df[[condition]]
      if (!is.null(controls)) {
        fmla <- as.formula(paste0('Surv(time_var, censor) ~ ', condition, ' + ', paste(controls, collapse = ' + ')))
      } else {
        fmla <- as.formula(paste0('Surv(time_var, censor) ~ ', condition))
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
            hazard_ratio = signif(exp(tidy_cox[1, 'estimate']), 4),
            ci_lower = signif(exp(tidy_confint[1, 1]), 4),
            ci_upper = signif(exp(tidy_confint[1, 2]), 4),
            pvalue = signif(tidy_cox[1, 'p.value'], 4)
          )
    }
  } else {
    # one hot encode the condition
    message('one hot encoding condition')
    df_ <- one_hot_encode_ovr(df, condition, binary = FALSE)
    vals <- unique(df_[, condition])
    columns_ovr <- colnames(df_)[colnames(df_) %in% paste0(condition, '_', vals)]

    # lapply the function over the columns_ovr and censors
    HR_list <- purrr::map(columns_ovr, ~ hazard_ratios_table(df_, .x, censors, controls = controls, censor_prefix = censor_prefix, time_prefix = time_prefix, ovr = FALSE))
  }

  HR_df <- do.call(rbind, HR_list)
  rownames(HR_df) <- gsub('censor_', '', rownames(HR_df))
  rownames(HR_df) <- gsub('\\.', '__', rownames(HR_df))
  return(HR_df)
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