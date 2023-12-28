################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Function to tabulate and summarize an input data - will return counts and also the stats testing the significance of the diff between groups

# Get the functions from github
library(devtools)
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/mgc_plotting_functions.R')


## EXAMPLE
# summarystatfile_seq <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2020/data/summarystats_IMGDEGIS_20200123.csv"
# sumstattab_seq <- summarize_table(intable = meta1[meta1[,"sequencing"] %in% "yes",], groupvar = "IMGDEGIS",
#                                   outfile = summarystatfile_seq, calc_stats = FALSE)



## A better summary function that will summarize table stats - first column should in general be some kind of ID
# Input: intable - the intable (first col is probably some kind of ID)
# summarystatfile - the file to write out too, this is a messy output, but is pretty good considering..
# will write out the file, and also output the summary list with each column as an entry in thelist

#### This is horrible coding, have to clean up at some point, but is working as intended now, so oh well.
summarize_table <- function(intable, groupvar = NULL, outfile, calc_stats = FALSE, calc_ci = FALSE) {

    ## Assign our infotab as the intable either as is, or with the grouping column excluded
    if (!is.null(groupvar)){
        groupcol <- intable[,groupvar,drop=FALSE]
        infotab <- intable[,!colnames(intable) %in% groupvar]
        subgroupvars <- sort(as.character(unique(groupcol[,groupvar])))
    } else {
        subgroupvars <- "all"
        infotab <- intable
    }

    ## For each column, we have to run analysis
    summarystatlist <- list()
    for (columnnum in seq_len(ncol(infotab))){

        ## For each subgroup, we have to tabulate our info
        ## If there is no goruping column, then we just run through once
        outstatlist <- list()
        for (groupvarnum in seq_len(length(subgroupvars))) {
            ## Grab the col
            if (!is.null(groupvar)){
                subtab <- infotab[groupcol[,groupvar] %in% subgroupvars[groupvarnum],columnnum,drop = FALSE]
            } else {
                subtab <- infotab[,columnnum,drop=FALSE]
            }
            statlabel <- colnames(subtab)

            ## This is where we actually do the tabulation
            ## if first column - thats the ID column, so just return the total and the unique
            ## if not, then see if column is character - then count, if number - then return statistical summary
            if (columnnum == 1){
                outstatlist[[groupvarnum]] <- cbind(classes = c("total", "unique_entries"), values = c(nrow(subtab), length(unique(subtab[,1]))))
            } else {
                if (class(subtab[,1]) == "character" | class(subtab[,1]) == "factor" | class(subtab[,1]) == "logical"){
                    ## failsafe for when the second group is all NAs
                    if (sum(is.na(subtab[,1])) == length(subtab[,1])) {
                        outstatlist[[groupvarnum]] <- data.frame(classes = "NAvalue", values = sum(is.na(subtab[,1])))
                    } else {
                        outstatlist[[groupvarnum]] <- data.frame(table(subtab[,1]))
                        colnames(outstatlist[[groupvarnum]]) <- c("classes", "values")
                    }
                }
                if (class(subtab[,1]) == "numeric" | class(subtab[,1]) == "integer") {
                    tempsum <- c(summary(subtab[,1]), sd = sd(subtab[,1], na.rm = TRUE))
                    outstatlist[[groupvarnum]] <- data.frame(ids = names(tempsum), values = as.numeric(unname(tempsum)))
                    colnames(outstatlist[[groupvarnum]]) <- c("classes", "values")
                }
            }
            names(outstatlist)[groupvarnum] <- subgroupvars[groupvarnum]
        }
        ## Merge our results back together
        outstat <- suppressWarnings(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "classes", all = TRUE, sort = FALSE), c(outstatlist)))
        outstat[is.na(outstat)] <- 0

        ## Changing the naming convention if we are grouping by a category
        if (is.null(groupvar)){
            colnames(outstat) <- c("classes", subgroupvars)
        } else {
            colnames(outstat) <- c(paste0(groupvar, "_cat"), subgroupvars)
        }

        ## Calculate the stats on our comparison (if applicable)
        if (calc_stats == TRUE) {
            if (columnnum == 1) {
                outstat <- outstat
            } else {
                if (class(infotab[,columnnum]) == "character" | class(infotab[,columnnum]) == "factor") {
                    if (length(subgroupvars) == 2) {
                        stattab <- outstat[,c(2:ncol(outstat))]
                        rownames(stattab) <- outstat[,1]

                        statval <- apply(stattab, 1, function(x) {
                            group1 <- c(rep(1,x[1]), rep(0, unname(colSums(stattab))[1] - x[1]))
                            group2 <- c(rep(1,x[2]), rep(0, unname(colSums(stattab))[2] - x[2]))
                            if(mean(group1) == mean(group2)) {
                                1
                            } else {
                                ## failsafe for case where you have 2 separate pops compared (all 1s vs all 0s)
                                ## Also need failsafe for when you have a group of size 1
                                
                                if (length(unique(group1)) == 1 & length(unique(group2)) == 1 | 
                                    length(group1) == 1 | length(group2) == 1) {
                                    pvalout <- 0
                                    if (calc_ci == TRUE) {
                                        ci_val <- "0 (0 - 0)"
                                        c(pvalout, ci_val)
                                    } else {
                                        pvalout
                                    }
                                } else {
                                    pvalout <- t.test(group1, group2)$p.value
                                    ## 20220228 - adding in a CI?
                                    if (calc_ci == TRUE) {
                                        grab_ci <- tryCatch(t.test(group1, group2)$conf.int,
                                                            error = function(e) c(0, 0))
                                        diff_val <- t.test(group1, group2)$estimate[1] / t.test(group1, group2)$estimate[2]
                                        ci_val <- paste0(diff_val, " (", paste(grab_ci, collapse = " - "), ")")
                                        c(pvalout, ci_val)
                                    } else {
                                        pvalout
                                    }
                                }
                            }
                        })
                        # outstat[,"statval"] <- statval
                        if (calc_ci == TRUE) {
                            outstat[,c("statval", "ci_val")] <- t(statval)
                        } else{
                            outstat[,"statval"] <- statval    
                        }
                        outstat <- outstat[order(factor(outstat[,1], levels = outstat[,1][order(as.character(outstat[,1]))])),]
                    } else {
                        ## NOT WRITTEN YET, WILL NEED TO WRITE AS APPLICABLE
                        # NEED TO DO A CHISQUARED TEST HERE!!!!
                        # aovtab <- intable[groupcol[,groupvar] %in% subgroupvars,c(statlabel, groupvar),drop = FALSE]
                        # aovtab[sapply(aovtab, is.character)] <- lapply(aovtab[sapply(aovtab, is.character)], as.factor)
                        # aovtab[sapply(aovtab, is.factor)] <- lapply(aovtab[sapply(aovtab, is.factor)],
                        #                                             function(x) `levels<-`(addNA(x), c(levels(x), "NA")))
                        # statval <- kruskal.test(as.formula(paste0(groupvar, "~", ".")), aovtab)$p.value
                        # statval <- chisq.test(groupcol[,1], infotab[,columnnum], simulate.p.value = TRUE)$p.value
                        statval <- tryCatch(chisq.test(groupcol[,1], infotab[,columnnum], simulate.p.value = TRUE)$p.value,
                                 error = function(e) NA)
                        outstat[,"statval"] <- statval
                        outstat <- outstat[order(factor(outstat[,1], levels = outstat[,1][order(as.character(outstat[,1]))])),]
                    }
                }
                if (class(infotab[,columnnum]) == "numeric" | class(infotab[,columnnum]) == "integer") {
                    ## if two, then t test
                    if (length(subgroupvars) == 2) {
                        group1 <- infotab[groupcol[,groupvar] %in% subgroupvars[1],columnnum,drop = FALSE]
                        group2 <- infotab[groupcol[,groupvar] %in% subgroupvars[2],columnnum,drop = FALSE]
                        if (length(group1[!is.na(group1)]) < 3 | length(group2[!is.na(group2)]) < 3) {
                            statval <- "toofew"
                            ci_val <- "toofew"
                        } else {
                            ## failsafe for case where you have 2 separate pops compared (all 1s vs all 0s)
                            if (nrow(unique(group1)) == 1 & nrow(unique(group2)) == 1) {
                                statval <- 0
                            } else {
                                ## Need a tryCatch here too for when you have data that is all 1s and NAs vs eachother (or other equiv scenario
                                statval <- tryCatch(t.test(group1, group2)$p.value,
                                                    error = function(e) NA)
                                # statval <- t.test(group1, group2)$p.value
                                ## 20220228 - adding in a CI?
                                grab_ci <- tryCatch(t.test(group1, group2)$conf.int,
                                                    error = function(e) c(0, 0))
                                diff_val <- t.test(group1, group2)$estimate[1] / t.test(group1, group2)$estimate[2]
                                ci_val <- paste0(diff_val, " (", paste(grab_ci, collapse = " - "), ")")
                            }
                        }
                        outstat[,"statval"] <- c(statval, rep("", nrow(outstat)-1))
                        if (calc_ci == TRUE) {outstat[,"ci_val"] <- c(ci_val, rep("", nrow(outstat)-1))}
                    } else {
                        ## NOT WRITTEN YET, WILL NEED TO WRITE AS APPLICABLE
                        # # if three or more, melt down into groups, and then anova test (melt so 2 cols, 1 with group 1 with stat)
                        # res.aov <- aov(var ~ group, data = melted data
                        # summary(res.aov)
                        aovtab <- intable[groupcol[,groupvar] %in% subgroupvars,c(statlabel, groupvar),drop = FALSE]
                        # an anova is finicky and needs to have the levels be real values (which isnt surprising) - so issues will happen with blanks and NAs, and there may be times where there are blanks AND NAs, so for this case, lets change them
                        aovtab[,groupvar][aovtab[,groupvar] == ""] <- "blankval"
                        aovtab[,groupvar][is.na(aovtab[,groupvar])] <- "NAval"
                        aovtab <- na.omit(aovtab)

                        ## Failsafe to get rid of NA values and turn into string which can then be factors
                        #`levels<-`(addNA(aovtab[,1]), c(levels(aovtab[,1]), "NA"))
                        aovtab[sapply(aovtab, is.character)] <- lapply(aovtab[sapply(aovtab, is.character)], as.factor)
                        # aovtab[sapply(aovtab, is.factor)] <- lapply(aovtab[sapply(aovtab, is.factor)],
                        #                                             function(x) `levels<-`(addNA(x), c(levels(x), "NA")))
                        if (nrow(aovtab) > 2) { ## Failsafe added 2021-9-14
                            # aov formula doesnt like dashes! So we have to change that: 2022-05-02
                            statlabel <- gsub("-", "_", statlabel)
                            colnames(aovtab)[1] <- gsub("-", "_", statlabel)
                            
                            aovform <- as.formula(paste0(statlabel, "~", "."))
                            aovout <- aov(aovform, aovtab)
                            statval <- summary(aovout)[[1]][["Pr(>F)"]][[1]]
                            outstat[,"statval"] <- statval
                            outstat <- outstat[order(factor(outstat[,1], levels = outstat[,1][order(as.character(outstat[,1]))])),]
                        } else {
                            outstat[,"statval"] <- NA
                        }

                    }
                }
            }
        }
        
        # outstat <- outstat[order(factor(outstat[,1], levels = outstat[,1][order(as.character(outstat[,1]))])),]
        
        ## Added a sorting function in here to help clean things up
        ## Write out the statlist and append to the table
        summarystatlist[[columnnum]] <- outstat
        names(summarystatlist)[columnnum] <- statlabel
        
        ## Adding a blank line to make the table a little easier to read?
        outstat <- rbind(outstat, rep(NA,nrow(outstat)))

        write.table(statlabel, outfile, sep = ",", append = columnnum != 1, row.names = FALSE, col.names = FALSE, quote = FALSE, na = "")
        suppressWarnings(write.table(outstat, outfile, sep = ",", append = TRUE, row.names = FALSE, col.names = !is.null(groupvar), quote = FALSE, na = ""))

    }
    return(summarystatlist)

}



## Play coding - turning my summary table output into a clean table using kable and a lot of formatting
## Ok - so we need to think about things that we want to do to make it easier
# 1 - move categories into the first column
# 2 - add percentages to the values
# 3 - round the stat values
# 4 - combining all into one big table

## NOT DOING the following
# picking which features to keep (just list all, should be easy to edit later)
# 
# sumstattable_input = sumstattab_seq
# addpercents = "vertical" # adds percents by category (as opposed to by feature)
# contsummary = c("mean", "sd") # what stats to pick for the continuous variable summary
# roundpvaldigits = 3
# c("median", "iqr")

clean_summarize_table_output <- function(sumstattable_input,
                                         addpercents = "vertical", contsummary = c("mean", "sd"), roundpvaldigits = 3) {
    ## Run over every entry in the list, and then clean as applicable:
    outtable_list <- list()
    for (stattablecatnum in seq_len(length(sumstattable_input))) {
        stattablecat_label <- names(sumstattable_input)[stattablecatnum]
        stattable_sel <- sumstattable_input[[stattablecatnum]]
        subgroups <- colnames(stattable_sel)[2:ncol(stattable_sel)][ 
            !grepl("statval|ci_val", colnames(stattable_sel)[2:ncol(stattable_sel)])]
        
        ## Ok, so the first one SHOULD be the ID number, and if it is, grab these values for later and move on
        if (identical(stattable_sel[,1], c("total", "unique_entries"))) {
            cohorttotals <- stattable_sel[stattable_sel[,1] == "total",2:ncol(stattable_sel)]
            cohorttotals_replacementtable <- data.frame(values = t(cohorttotals),
                                                        original = names(cohorttotals),
                                                        replacement = paste0(names(cohorttotals), " (N = ", cohorttotals, ")"))
            next
        }
        
        ## Next - if we have a continuous variable, then we will have the same suumary stats, so we can detect by that, and then adjust accordingly
        if (sum(stattable_sel[,1] %in% c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "sd")) >= 7) { ## we have a cont variable
            # if(contsummary[1] == "mean") {var1 <- stattable_sel[stattable_sel[,1] == "Mean",2:(ncol(stattable_sel)-1)]}
            # if(contsummary[1] == "median") {var1 <- stattable_sel[stattable_sel[,1] == "Median",2:(ncol(stattable_sel)-1)]}
            # if(contsummary[2] == "sd") {var2 <- stattable_sel[stattable_sel[,1] == "sd",2:(ncol(stattable_sel)-1)]}
            # if(contsummary[2] == "iqr") {var2 <- stattable_sel[stattable_sel[,1] == "1st Qu.",2:(ncol(stattable_sel)-1)]}
            if(contsummary[1] == "mean") {var1 <- stattable_sel[stattable_sel[,1] == "Mean",subgroups]}
            if(contsummary[1] == "median") {var1 <- stattable_sel[stattable_sel[,1] == "Median",subgroups]}
            if(contsummary[2] == "sd") {var2 <- round(stattable_sel[stattable_sel[,1] == "sd",subgroups],2)}
            if(contsummary[2] == "iqr") {var2 <- paste0(
                round(stattable_sel[stattable_sel[,1] == "1st Qu.",subgroups], 2), "-", 
                round(stattable_sel[stattable_sel[,1] == "3rd Qu.",subgroups],2))
                                                        }
            outvar <- paste0(round(var1,2), " (", var2, ")")
            names(outvar) <- names(var1)
            
            outstat <- suppressWarnings(round(as.numeric(stattable_sel[1,"statval"]), roundpvaldigits))
            if (length(outstat) > 0){
                if (outstat %in% 0) {outstat <- "<0.001"}
            } else{
                outstat <- NULL
            }

            outtable <- cbind(
                # category = paste0(stattablecat_label, ", ", contsummary[1], " (±", contsummary[2], ")"),
                category = stattablecat_label,
                feature = paste0(contsummary[1], " (±", contsummary[2], ")"),
                t(data.frame(outvar)),
                statval = outstat)
            if (ncol(outtable) == 3) {colnames(outtable) <- c("category", "feature", "all")}
        } else { ## THIS SHOULD BE CATEGORICAL, I CANT THINK OF ANYTHING ELSE THIS COULD BE RIGHT NOW...
            outtable <- cbind(category = stattablecat_label, stattable_sel)
            
            ## If we have a multigroup summary with stats do this:
            if (sum(grepl("statval", colnames(stattable_sel))) > 0) {
                colnames(outtable) <- c("category", "feature", colnames(stattable_sel)[2:(ncol(stattable_sel)-1)], "statval")
                outtable[,"statval"] <- round(as.numeric(outtable[,"statval"]), roundpvaldigits)
                outtable[,"statval"] <- ifelse(outtable[,"statval"] == 0, "<0.001", outtable[,"statval"])
            } else {
                ## Lastly - if we have a single group summary and no stats, then do this:
                if (ncol(stattable_sel) == 2){
                    outtable <- cbind(category = stattablecat_label, stattable_sel)
                    colnames(outtable) <- c("category", "feature", "all")
                }  else { ## If we have multigroup with NO stats do this:
                    outtable <- cbind(category = stattablecat_label, stattable_sel)
                    colnames(outtable) <- c("category", "feature", colnames(stattable_sel)[2:ncol(stattable_sel)])
                }
            }

            ## Now add some things
            if (!is.null(addpercents)) {
                outtable[,"feature"] <- paste0(outtable[,"feature"], ", N (%)")
                ## percents by row or column?
                percentorderparam <- ifelse(addpercents == "vertical", 2, 1)
                # outtable[,colnames(stattable_sel)[2:(ncol(stattable_sel)-1)]] <- apply(
                #     outtable[,colnames(stattable_sel)[2:(ncol(stattable_sel)-1)]], percentorderparam, function(x) {
                #         percentvals <- round(x/sum(x),3)*100
                #         outvals <- paste0(x, " (", percentvals, ")")
                #         outvals
                #     }) ## TOSH CHANGED - the minus one correction needs to be changed, just want to grab all groups - all for no groups
                tempout <- apply(outtable[,subgroups,drop=FALSE], percentorderparam, function(x) {
                        percentvals <- round(x/sum(x),3)*100
                        outvals <- paste0(x, " (", percentvals, ")")
                        outvals
                    })
                if (percentorderparam == 1) {tempout <- t(tempout)}
                outtable[,subgroups] <- tempout
            } else{
                outtable[,"feature"] <- paste0(outtable[,"feature"], ", N")
            }
        }
        outtable_list[[stattablecatnum]] <- outtable
    }
    combotable <- do.call(rbind, outtable_list)
    rownames(combotable) <- NULL
    colnames(combotable)[3:(ncol(sumstattable_input[[1]])+1)] <- cohorttotals_replacementtable$replacement[match(colnames(combotable)[3:(ncol(sumstattable_input[[1]])+1)], cohorttotals_replacementtable$original)]
    return(combotable)
}




## A better summary function that will summarize table stats - first column should in general be some kind of ID
# Input: intable - the intable (first col is probably some kind of ID)
# summarystatfile - the file to write out too, this is a messy output, but is pretty good considering..
# will write out the file, and also output the summary list with each column as an entry in thelist

#### This is horrible coding, have to clean up at some point, but is working as intended now, so oh well.
summarize_and_correlate_table <- function(intable, groupvar) {
    
    ## Assign our infotab as the intable either as is, or with the grouping column excluded
    groupcol <- intable[,groupvar,drop=FALSE]
    infotab <- intable[,!colnames(intable) %in% groupvar]
    
    ## For each column, we have to run analysis
    outstatlist <- list()
    for (columnnum in seq_len(ncol(infotab))){
    # for (columnnum in seq_len(11)){
        subtab <- infotab[,columnnum,drop=FALSE]
        statlabel <- colnames(subtab)

        ## If first column - then assume IDs, and just count them up
        if (columnnum == 1){
            outstatlist[[columnnum]] <- data.frame(rval = nrow(subtab), pval = length(unique(subtab[,1])), row.names = statlabel)
        } else {
            if (class(subtab[,1]) == "character" | class(subtab[,1]) == "factor"){
                ## NOT WRITTEN YET, WILL NEED TO WRITE AS APPLICABLE
                # # if three or more, melt down into groups, and then anova test (melt so 2 cols, 1 with group 1 with stat)
                # res.aov <- aov(var ~ group, data = melted data
                # summary(res.aov)
                aovtab <- intable[,c(statlabel, groupvar),drop = FALSE]
                # an anova is finicky and needs to have the levels be real values (which isnt surprising) - so issues will happen with blanks and NAs, and there may be times where there are blanks AND NAs, so for this case, lets change them
                if (length(aovtab[,statlabel][aovtab[,statlabel] == ""]) > 0 ) {
                    levels(aovtab[,statlabel]) <- c(levels(aovtab[,statlabel]), "blankval")
                    aovtab[,statlabel][aovtab[,statlabel] == ""] <- "blankval"
                }
                if (length(aovtab[,statlabel][is.na(aovtab[,statlabel])]) > 0 ) {
                    levels(aovtab[,statlabel]) <- c(levels(aovtab[,statlabel]), "NAval")
                    aovtab[,statlabel][is.na(aovtab[,statlabel])] <- "NAval"
                }
                aovtab <- na.omit(aovtab)
                ## In case these values are characters, need to be factors for anova
                aovtab[sapply(aovtab, is.character)] <- lapply(aovtab[sapply(aovtab, is.character)], as.factor)

                aovform <- as.formula(paste0(groupvar, "~", statlabel))
                aovout <- aov(aovform, aovtab)
                statval <- summary(aovout)[[1]][["Pr(>F)"]][[1]]
                outstatlist[[columnnum]] <- data.frame(rval = 0, pval = statval, row.names = statlabel)
            }
            if (class(subtab[,1]) == "numeric" | class(subtab[,1]) == "integer") {
                suppressWarnings(corout <- cor.test(groupcol[,1], subtab[,1]))
                outstatlist[[columnnum]] <- data.frame(rval = corout$estimate, pval = corout$p.value, row.names = statlabel)
            }
        }
        names(outstatlist)[columnnum] <- statlabel
    }
    outstattable <- do.call(rbind, outstatlist)
    return(outstattable)
    
}












# if (columnnum == 1){
#     outstatlist[[groupvarnum]] <- cbind(classes = c("total", "unique_entries"), values = c(nrow(subtab), length(unique(subtab[,1]))))
# } else {
#     if (class(subtab[,1]) == "character" | class(subtab[,1]) == "factor"){
#         outstatlist[[groupvarnum]] <- data.frame(table(subtab[,1]))
#         colnames(outstatlist[[groupvarnum]]) <- c("classes", "values")
#     }
#     if (class(subtab[,1]) == "numeric" | class(subtab[,1]) == "integer") {
#         tempsum <- summary(subtab[,1])
#         outstatlist[[groupvarnum]] <- data.frame(ids = names(tempsum), values = as.numeric(unname(tempsum)))
#         colnames(outstatlist[[groupvarnum]]) <- c("classes", "values")
#     }
# }




## Generate Some prestats for this analysis
### IF YOU WANT TO FIX THIS TO OUTPUT GGPLOTS - THEN USE THIS
### https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r

# THE FIRST COL WILL BE SAMPLE IDS TO MATCH THE FUNCTION ABOVE

# intable - the table with all the information to match from, and rownames being the samples
# groupcat - the label to split up by groups
# controlcat - the labels that are the categories to summarize
# discretevar - the sublabels of controlcat that are discrete
# contvar - the sublabels of controlcat that are continuous
summarize_table_figures <- function(intable, groupvar = NULL, outfilepath) {

    dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)

    ## Assign our infotab as the intable either as is, or with the grouping column excluded
    if (!is.null(groupvar)){
        # REMOVE THE SAMPLEID COLUMN THAT IS THE FIRST COLUMN
        infotab <- intable[,-1]
        ## Maybe I need to turn any NA in the group into an "NAgroup" so they still get counted as theyre own thing
        infotab[is.na(infotab[,groupvar]),groupvar] <- "NAgroup"
        subgroupvars <- as.character(unique(infotab[,groupvar]))
    } else {
        subgroupvars <- "ALLSAMPLES"
        groupvar <- "ALLSAMPLES"
        # REMOVE THE SAMPLEID COLUMN THAT IS THE FIRST COLUMN
        infotab <- cbind(ALLSAMPLES = "ALLSAMPLES", intable[,-1])
    }

    ## For each column, we have to run analysis
    for (columnnum in seq_len(ncol(infotab)-1)){
        catlabel <- colnames(infotab)[columnnum+1]
        subtable <- infotab[,c(groupvar, catlabel)]

        # If discrete - lets do stack bar charts with % and number for each group
        if (class(subtable[,2]) == "character" | class(subtable[,2]) == "factor") {
            meltedtab  <- melt(table(subtable))
            meltedtab[,2] <- as.character(meltedtab[,2])

            # stack bar chart with N
            pout <- ggplot(meltedtab, aes(fill = meltedtab[,catlabel], y = meltedtab[,"value"], x = meltedtab[,groupvar]))
            pout <- pout + geom_bar(position="stack", stat = "identity")
            pout <- pout + labs(title = catlabel, x = "category", y = "value", fill = catlabel)
            pout <- pout + theme_pubr(base_size = 10, x.text.angle = 60)
            pdf(paste0(outfilepath, catlabel, "_stack_barN.pdf"))
            print(pout)
            junk <- dev.off()

            # stack bar chart with %
            pout2 <- ggplot(meltedtab, aes(fill = meltedtab[,catlabel], y = meltedtab[,"value"], x = meltedtab[,groupvar]))
            pout2 <- pout2 + geom_bar(position="fill", stat = "identity")
            pout2 <- pout2 + labs(title = catlabel, x = "category", y = "value", fill = catlabel)
            pout2 <- pout2 + theme_pubr(base_size = 10, x.text.angle = 60)
            pdf(paste0(outfilepath, catlabel, "_stack_barperc.pdf"))
            print(pout2)
            junk <- dev.off()

        }

        # If continuous - lets do histograms
        if (class(subtable[,2]) == "numeric" | class(subtable[,2]) == "integer") {

            ## Histograms plots
            labsparam <- list(title = paste0(catlabel, " histogram"), x = catlabel, y = "count")
            histout <- plot_histogram(data = subtable[,2,drop=FALSE], groupvar = subtable[,groupvar,drop=FALSE],
                                      fitcurve = TRUE, labsparam, limitx = NULL, limity = NULL,
                                      # binparam=(max(subtable[,2], na.rm = TRUE) - min(subtable[,2], na.rm = TRUE)))
                                      binparam=min(length(unique(subtable[,2])), 20))
            ## Boxplots
            labsparam <- list(title = paste0(catlabel, " boxplot"), x = catlabel, y = "value")
            plotstatparam <- ifelse(length(unique(subtable[,groupvar])) > 1, "intra", FALSE)
            bpintable <- melt(subtable, id.vars = groupvar)
            bpout <- boxplot_plotter(boxplottable = bpintable, xsplit = "category", labsparam = labsparam,
                                     testtypeparam = "t.test", colorparam = NULL, secondaxis=NULL,
                                     plotstats = plotstatparam, comparisonparam = unique(subtable[,groupvar]))

            ## Write out plots
            pdf(paste0(outfilepath, catlabel, "_hist.pdf"))
            print(histout)
            junk <- dev.off()
            pdf(paste0(outfilepath, catlabel, "_bp.pdf"))
            print(bpout)
            junk <- dev.off()

        }
    }
}





