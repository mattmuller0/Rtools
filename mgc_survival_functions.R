################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Survival Analysis Functions

## Load in Libraries
packagelist = c("ggplot2", "reshape2", "DESeq2", "grid", "gridExtra", "scales", "ggrepel", "tools", "Hmisc", "survminer", "plyr")
junk <- lapply(packagelist, function(xxx) suppressMessages(
    require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))



## Input to functions needs to be survival data:
# THERE REALLY ISNT ANY WAY TO ADD COVARIATES FOR THIS - YOU SHOULD JUST CREATE YOUR COHORTS SEPARATELY, AND ADD IT TO THE COHORT COL 
## EX) (COH1_MALE, COH2_MALE, COH1_FEMALE, COH2_FEMALE, ETC.)
# table with rownames - then (1) Time to event, (2) event status, (3) the cohort that it is in.
# time breakparam will manually split by the value you give, otherwise - will break into increments of 10
create_survival_plot <- function(survivaldata, timebreakparam = NULL, ylimitparam = NULL, confintervalparam = FALSE) {
    
    ## Read in the table and grab all the original labels (need for formula later)
    ## Rename with specific names - necessary for all of the formulas you need to do.
    modeltable <- survivaldata
    colnames(modeltable) <- c("Time_to_event", "event_status", "Cohort")

    ## Set a proper time break param
    if (is.null(timebreakparam)) {
        maxtimeparam <- max(modeltable[,"Time_to_event"])
        timebreakparam <- ceiling(max(modeltable[,"Time_to_event"])/10)
    } else {
        maxtimeparam <- round_any(max(modeltable[,"Time_to_event"]), 100, f = ceiling)
    }
    
    ## Create the sruvival object
    survobject <- surv_fit(Surv(Time_to_event, event_status) ~ Cohort, data = modeltable)
    
    ## Save out the survival p value for output
    outsurvpvalue <- surv_pvalue(survobject)
    ## Create the table of data
    res <- summary(survobject, times = seq(0,timebreakparam*ceiling(maxtimeparam/timebreakparam), timebreakparam))
    outsurvtable <- as.data.frame(res[c("strata", "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper", "cumhaz", "std.chaz")])
    
    ## KM curves with a log-rank p value do NOT have any intrinsic effect sizes, so we have to resort to another analysis (coxph) to get that info
    # TOSH - DO YOU NEED TO ADD A REF PARAM HERE???
    # Suppressing warnings here due to inf coefficient warning 2021/02/02
    suppressWarnings(outcoxphobject <- coxph(Surv(Time_to_event, event_status) ~ Cohort, data = modeltable))

    ## SURVIVAL CURVES CANT WORK WITH CONTINUOUS VARIABLES, BUT COXPH CAN
    
    ## Added a ylimparam for smaller changes that I need to expand to plot - this requires a param for pval.coord change as well
    if (!is.null(ylimitparam)) {
        pvalcordparam <- c(max(modeltable[,1]) * 0.01, min(ylimitparam) + diff(range(ylimitparam))/10 * 3)
        pvalmethodcordparam <- c(max(modeltable[,1]) * 0.01, min(ylimitparam) + diff(range(ylimitparam))/10 * 1)
    } else {
        pvalcordparam <- pvalmethodcordparam <- NULL
    }
   
    
    outsurvplot <- ggsurvplot(
                              fit = survobject,
                              data = modeltable,
                              risk.table = TRUE,
                              conf.int = confintervalparam,
                              conf.int.style = "ribbon",
                              pval = TRUE, pval.method = TRUE,
                              break.time.by = timebreakparam,
                              ggtheme = theme_minimal(base_size = 10),
                              risk.table.y.text.col = TRUE,
                              risk.table.y.text = FALSE,
                              ylim = ylimitparam, pval.coord = pvalcordparam, pval.method.coord = pvalmethodcordparam)

    return(list(outsurvtable = outsurvtable, outsurvplot = outsurvplot, 
                outsurvpvalue = outsurvpvalue, outsurvobject = survobject, 
                outcoxphobject = outcoxphobject))
    
}


## Input to functions needs to be survival data:
# table with rownames - then (1) Time to event, (2) event status, (3) the cohort that it is in.
# optional parameter for the covariates - WHICH APPARENTLY IS JUST USELESS, YOU JUST ADD ADDITIONAL COLUMNS TO THE SURVIVAL INDATA
coxph_analysis <- function(survivaldata) {

    ## Rename with specific names - necessary for all of the formulas you need to do.
    modeltable <- survivaldata
    colnames(modeltable)[1:2] <- c("Time_to_event", "event_status")
    ## Failsafe to make column names formula safe:
    colnames(modeltable) <- gsub("-", "\\.", colnames(modeltable))
    
    coxph_object <- coxph(as.formula(paste0("Surv(Time_to_event, event_status) ~ ", paste(colnames(modeltable)[3:(ncol(modeltable))], collapse = " + "))), 
                          data =  modeltable)
    coxph_outtable <- merge(summary(coxph_object)[[7]], summary(coxph_object)[[8]], by = c("row.names", "exp(coef)"))
    
    # I need a failsafe in here: 09/28/2022
    ## Heres a thought - a failsafe check that asks if any covar has 0 events, that combine that level with the next biggest level
    # if (sum(unlist(coxph_outtable[,grepl("\\.95", colnames(coxph_outtable))] > 100)) > 0) {  # If there are any "Inf" or massive CI vals, then remove that label
    #     
    # }
    
    # covarlabels <- colnames(survivaldata)[3:ncol(survivaldata)]
    # for (covarnum in seq_len(length(covarlabels))) {coxph_outtable[,1] <- gsub(paste0("covar", covarnum), paste0(covarlabels[covarnum], "__"), coxph_outtable[,1])}
    coxph_outplot <- ggforest(coxph_object, data = modeltable, fontsize = 1)

    
    ## Read in the table and grab all the original labels (need for formula later)
    # intable <- survivaldata
    # inlabels <- colnames(survivaldata)
    # intime <- inlabels[1]
    # inevent <- inlabels[2]
    # incohort <- inlabels[3]
    # colnames(modeltable) <- c("Time_to_event", "event_status", "Cohort")
    # colnames(modeltable)[3:(ncol(modeltable))] <- paste0("covar", seq(ncol(modeltable)-2))
    
    # coxph_object <- coxph(Surv(Time_to_event, event_status) ~ Cohort, data =  modeltable)
    # coxph_outtable <- merge(summary(coxph_object)[[7]], summary(coxph_object)[[8]], by = c("row.names", "exp(coef)"))
    # coxph_outplot <- ggforest(coxph_object, data = modeltable, fontsize = 1.3)
    # 
    # ## If there is covar data, read in the covariate data and merge it with out other table
    # if (!is.null(covardata)) {
    #     modeltable <- merge(modeltable, covardata, by = "row.names")
    #     rownames(modeltable) <- modeltable[,1]
    #     modeltable <- modeltable[,!grepl("Row.names", colnames(modeltable))]
    #     colnames(modeltable)[4:(3+ncol(covardata))] <- paste0("covar", seq(ncol(covardata)))
    #     
    #     coxph_object <- coxph(as.formula(paste0("Surv(Time_to_event, event_status) ~ Cohort + ", 
    #                                                paste(colnames(modeltable)[4:(3+ncol(covardata))], collapse = " + "))), 
    #                              data =  modeltable)
    #     coxph_outtable <- merge(summary(coxph_object)[[7]], summary(coxph_object)[[8]], by = c("row.names", "exp(coef)"))
    #     ## Put the original labels back
    #     covarlabels <- colnames(covardata)
    #     for (covarnum in seq_len(length(covarlabels))) {coxph_outtable[,1] <- gsub(paste0("covar", covarnum), covarlabels[covarnum], coxph_outtable[,1])}
    #     
    #     coxph_outplot <- ggforest(coxph_object, data = modeltable, fontsize = 1)
    # }
    
    return(list(coxph_outtable = coxph_outtable, coxph_outplot = coxph_outplot))

}

