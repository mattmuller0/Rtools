################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## ISCHEMIA Subtype Analysis

## Load in Libraries
packagelist = c("NMF", "Hmisc", "tools", "caret", "blacksheepr", "mlbench", "caret", "pROC", "WGCNA", "multiROC", "scattermore", "dendextend", "devtools")
junk <- lapply(packagelist, function(xxx) suppressMessages(
    require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/mgc_plotting_functions.R')
# source("/Users/tosh/Desktop/Ruggles_Lab/code/process_metadata_functions.R")
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/mgc_survival_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/mgc_WGCNA_functions.R')

# --------------------------------- Add on and edit the biorep table --------------------------------- 

addon_and_edit_biorep_table <- function(bioreptable, addclustertypes = TRUE, ...) {
    
    
    addontable1 <- data.frame(PATNUM = SOIall, CUSTOM_maxfollow = apply(bioreptable[SOIall,grepl("^T_", colnames(bioreptable))], 1, function(x) max(x)/365.25))
    addontable2 <- data.frame(PATNUM = SOIall, CUSTOM_BMI_obese = ifelse(bioreptable[,"BMI"] >= 30, "Obese", "NotObese"))
    addontable3 <- data.frame(PATNUM = SOIall, CUSTOM_SMOKSTAT_currentsmoker = ifelse(bioreptable[,"SMOKSTAT"] %in% "Current Smoker", "CurrentSmoker", "NotCurrentSmoker"))
    addontable4 <- data.frame(PATNUM = SOIall, CUSTOM_PCIorCABG = ifelse(!is.na(bioreptable[,"FIRSTRV"]), "PCIorCABG", "NoPCIorCABG"))
    # addontable5 <- data.frame(PATNUM = SOI, CUSTOM_PAD_NOCCEREB = bioreptable[,"PRCERPAD"] - bioreptable[,"PRICEREB"])
    addontable6 <- data.frame(PATNUM = SOIall, CUSTOM_RACE_AGGR = ifelse(bioreptable[,"RACE"] %in% c("White", "Black or African American", "Asian"), bioreptable[,"RACE"], "Other or multiple ethnic groups"))
    addontable7 <- data.frame(PATNUM = SOIall, CUSTOM_LVEFCONT_LVEFund55 = ifelse(bioreptable[,"LVEFCONT"] < 55, "LVEFund55", "NoLVEFund55"))
    addontable8 <- data.frame(PATNUM = SOIall, CUSTOM_EGFR_EGFRund60 = ifelse(bioreptable[,"EGFR"] < 60, "EGFRund60", "NoEGFRund60"))
    addontable9 <- data.frame(PATNUM = SOIall, CUSTOM_LDLC_LDLCund70 = ifelse(bioreptable[,"LDLC"] < 70, "LDLCund70", "NoLDLCund70"))
    addontable9_2 <- data.frame(PATNUM = SOIall, CUSTOM_LDLC_LDLCund100 = ifelse(bioreptable[,"LDLC"] < 100, "LDLCund100", "NoLDLCund100"))
    addontable10 <- data.frame(PATNUM = SOIall, CUSTOM_RVBPSYS_RVBPSYSund140 = ifelse(bioreptable[,"RVBPSYS"] < 140, "RVBPSYSund140", "NoRVBPSYSund140"))
    addontable11 <- data.frame(PATNUM = SOIall, CUSTOM_LDLCund70andMESTATIN = ifelse(bioreptable[,"LDLC"] < 70 & bioreptable[,"MESTATIN"] == "Yes", "LDLCund70andMESTATIN", "NoLDLCund70andMESTATIN"))
    addontable12 <- data.frame(PATNUM = SOIall, CUSTOM_CKD_TRT = ifelse(bioreptable[,"CKD"] %in% "Yes", bioreptable[,"TRT"], NA))
    addontable13 <- data.frame(PATNUM = SOIall, CUSTOM_ACEIandARB = ifelse(bioreptable[,"MEAHACE"] %in% "Yes" | bioreptable[,"MEAHARB"] %in% "Yes", "ACEIorARB", "NoACEIorARB"))
    addontable14 <- data.frame(PATNUM = SOIall, CUSTOM_DEGRISCH_NONEMILD = ifelse(bioreptable[,"DEGRISCH"] %in% c("None", "Mild"), "NoneMild", bioreptable[,"DEGRISCH"]))
    addontable15 <- data.frame(PATNUM = SOIall, CUSTOM_IMGDEGIS_NONEMILD = ifelse(bioreptable[,"IMGDEGIS"] %in% c("None", "Mild"), "NoneMild", bioreptable[,"IMGDEGIS"]))
    
    # addontable16 <- data.frame(PATNUM = SOIall, DUKESCORE_NUMERIC = ifelse(bioreptable[,"DUKESCORE"] %in% c("Mild", "None"), "NoneMild", bioreptable[,"IMGDEGIS"]))
    # 
    # 
    # tt1 <- align_metatable(intable = bioreptable[bioreptable[,"PATNUM"] %in% SOIall, "DUKESCORE", drop=FALSE], list(
    #     "DUKESCORE" = list("DUKESCORE_7" = "Left Main >=50%", 
    #                        "DUKESCORE_6" = "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD", 
    #                        "DUKESCORE_5" = "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
    #                        "DUKESCORE_4" = "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
    #                        "DUKESCORE_3" = "1 Vessel with at least Moderate (>=50%) Plaque")))
    
    
    # Adding on some cleaned columns of heart disease variables
    IMGDEGIS_addontable <- bioreptable[,c("PATNUM", "IMGDEGIS"), drop=FALSE]
    IMGDEGIS_addontable[,"IMGDEGIS_01_2_3"] <- ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% c("None", "Mild"), "NoneMild",
                                                      ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% "Moderate", "Moderate",
                                                             ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% "Severe", "Severe", NA)))
    IMGDEGIS_addontable[,"IMGDEGIS_01_3"] <- ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% c("None", "Mild"), "NoneMild",
                                                    ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% "Severe", "Severe", NA))
    
    CTNDV50_addontable <- bioreptable[,c("PATNUM", "CTNDV50"), drop=FALSE]
    CTNDV50_addontable[,"CTNDV50_123_clean"] <- ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "1", "1",
                                                       ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "2", "2",
                                                              ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "3", "3", NA)))
    CTNDV50_addontable[,"CTNDV50_13_clean"] <- ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "1", "1",
                                                      ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "3", "3", NA))
    
    CTNDV70_addontable <- bioreptable[,c("PATNUM", "CTNDV70"), drop=FALSE]
    CTNDV70_addontable[,"CTNDV70_0123_clean"] <- ifelse(CTNDV70_addontable[,"CTNDV70"] %in% "0", "0",
                                                ifelse(CTNDV70_addontable[,"CTNDV70"] %in% "1", "1",
                                                ifelse(CTNDV70_addontable[,"CTNDV70"] %in% "2", "2",
                                                ifelse(CTNDV70_addontable[,"CTNDV70"] %in% "3", "3", NA))))
    CTNDV70_addontable[,"CTNDV70_13_clean"] <- ifelse(CTNDV70_addontable[,"CTNDV70"] %in% "1", "1",
                                               ifelse(CTNDV70_addontable[,"CTNDV70"] %in% "3", "3", NA))
    CTNDV70_addontable[,"CTNDV70_combo0noneval"] <- ifelse(CTNDV70_addontable[,"CTNDV70"] %in% c("Non-evaluable", "0"), "0orNoneval", 
                                                    ifelse(CTNDV70_addontable[,"CTNDV70"] %in% c(""), NA, CTNDV70_addontable[,"CTNDV70"]))
    
    
    ## Also adding on the clustermembership tables to this:
    # rna_clustermembership_table, meth_clustermembership_table
    print(addclustertypes)
    if (addclustertypes) {
        rna_clustermembership_addon <- data.frame(PATNUM = SOIall, rna_clustermembership_table[SOIall, c("rna_4cluster_w3AB", "rna_4cluster")])
        meth_clustermembership_addon <- data.frame(PATNUM = SOIall, meth_clustermembership_table[SOIall, c("meth_3cluster", "meth_4cluster")])
    } else {
        rna_clustermembership_addon <- meth_clustermembership_addon <- NULL
    }

    
    bioreptable_waddons <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE), 
                                  Filter(Negate(is.null), list(  # remove the null entries0
        bioreptable,
        addontable1, addontable2, addontable3, addontable4, # addontable5,
        addontable6, addontable7, addontable8, addontable9, addontable9_2, addontable10, 
        addontable11, addontable12, addontable13, addontable14, addontable15,
        biomarkertable_sel,
        IMGDEGIS_addontable[,c("PATNUM", "IMGDEGIS_01_2_3", "IMGDEGIS_01_3")],
        CTNDV50_addontable[,c("PATNUM", "CTNDV50_123_clean", "CTNDV50_13_clean")],
        CTNDV70_addontable[,c("PATNUM", "CTNDV70_0123_clean", "CTNDV70_13_clean", "CTNDV70_combo0noneval")],
        rna_clustermembership_addon, meth_clustermembership_addon
    )))
    rownames(bioreptable_waddons) <- bioreptable_waddons[,"PATNUM"]
    
    ## Change the event columns to characters for proper tabulations:
    cols_to_characterify <- c(colnames(bioreptable_waddons)[grepl("^C_", colnames(bioreptable_waddons))], "PRICEREB")
    bioreptable_waddons[,cols_to_characterify] <- apply(bioreptable_waddons[,cols_to_characterify], 2, function(x) as.character(x))
    
    # Make all NA values set to NotAvailable
    # tt1 <- bioreptable_waddons
    
    return(bioreptable_waddons)
}


# --------------------------------- Factorize our biorep table with preset levels --------------------------------- 

# Create a factorized table for whenever that particular things matters (with presets ready to go)
factorize_metatable <- function(intable) {
    outtable <- intable
    COInames <- colnames(outtable)
    
    if ("IMGDEGIS" %in% COInames) {outtable[,"IMGDEGIS"] <- factor(outtable[,"IMGDEGIS"], levels = c("None", "Mild", "Moderate", "Severe"))}
    if ("IMGDEGIS_01_2_3" %in% COInames) {outtable[,"IMGDEGIS_01_2_3"] <- factor(outtable[,"IMGDEGIS_01_2_3"], levels = c("NoneMild", "Moderate", "Severe"))}
    if ("IMGDEGIS_01_3" %in% COInames) {outtable[,"IMGDEGIS_01_3"] <- factor(outtable[,"IMGDEGIS_01_3"], levels = c("NoneMild", "Severe"))}
    
    if ("CTNDV50" %in% COInames) {outtable[,"CTNDV50"] <- factor(outtable[,"CTNDV50"], levels = c("1", "2", "3", "Non-evaluable"))}
    if ("CTNDV50_123_clean" %in% COInames) {outtable[,"CTNDV50_123_clean"] <- factor(outtable[,"CTNDV50_123_clean"], levels = c("1", "2", "3"))}
    if ("CTNDV50_13_clean" %in% COInames) {outtable[,"CTNDV50_13_clean"] <- factor(outtable[,"CTNDV50_13_clean"], levels = c("1", "3"))}
    
    if ("CTNDV70" %in% COInames) {outtable[,"CTNDV70"] <- factor(outtable[,"CTNDV70"], levels = c("Non-evaluable", "0", "1", "2", "3"))}
    if ("CTNDV70_combo0noneval" %in% COInames) {outtable[,"CTNDV70_combo0noneval"] <- factor(outtable[,"CTNDV70_combo0noneval"], levels = c("0orNoneval", "1", "2", "3"))}
    
    if ("CKD" %in% COInames) {outtable[,"CKD"] <- factor(outtable[,"CKD"], levels = c("No", "Yes"))}
    if ("highlowage" %in% COInames) {outtable[,"highlowage"] <- factor(outtable[,"highlowage"], levels = c("AGE_RAND_low", "AGE_RAND_highQ3"))}
    
    if ("rna_4cluster_w3AB" %in% COInames) {outtable[,"rna_4cluster_w3AB"] <- factor(outtable[,"rna_4cluster_w3AB"],
        levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3A", "nmf_cluster_3B", "nmf_cluster_4"))}
    if ("rna_4cluster" %in% COInames) {outtable[,"rna_4cluster"] <- factor(outtable[,"rna_4cluster"],
        levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4"))}
    if ("meth_3cluster" %in% COInames) {outtable[,"meth_3cluster"] <- factor(outtable[,"meth_3cluster"],
        levels = c("meth_3cluster_1", "meth_3cluster_2", "meth_3cluster_3"))}
    if ("meth_4cluster" %in% COInames) {outtable[,"meth_4cluster"] <- factor(outtable[,"meth_4cluster"],
        levels = c("meth_4cluster_1", "meth_4cluster_2", "meth_4cluster_3", "meth_4cluster_4"))}
    
    if ("CAGE_RND" %in% COInames) {outtable[,"CAGE_RND"] <- factor(outtable[,"CAGE_RND"], 
        levels = c("<= 34", "35 - 44", "45 - 54", "55 - 64", "65 - 74", ">= 75"))}
    if ("SEX" %in% COInames) {outtable[,"SEX"] <- factor(outtable[,"SEX"], levels = c("Female", "Male"))}
    if ("RACE" %in% COInames) {outtable[,"RACE"] <- factor(outtable[,"RACE"], 
        levels = c("White", "American Indian or Alaska Native", "Black or African American", "Asian", "Native Hawaiian or Other Pacific Islander",
                   "Multiple Races", ""))}
    
    return(outtable)
}

# --------------------------------- custom heatmap color annotation list from the color guide and COI --------------------------------- 

# Create heatmap annotation custom list from the colorguide and selected COI
custom_annotation_list_from_colorguide <- function(COI, colorguide) {
    COIlist <- list()
    for (COInum in seq_len(length(COI))) {
        COIsel <- COI[COInum]
        COIlist[COInum] <- list(colorguide[colorguide[,"category"] %in% COIsel,"color"])
        names(COIlist[[COInum]]) <- colorguide[colorguide[,"category"] %in% COIsel,"feature"]
        names(COIlist)[COInum] <- COIsel
    }
    return(COIlist)
}


# --------------------------------- ischemia disease metric plot --------------------------------- 
# Second viz - the overlap of the ischemia/CT variables - I really think that is a useful viz
ischemia_disease_metrics_plot <- function(bioreptable_waddons, colorguide, plotmetrics = c("IMGDEGIS_01_2_3", "CTNDV70_combo0noneval")) {
    
    # create a factorized metatable for our variables of interest, factors preserve order
    labeltable_temp <- factorize_metatable(bioreptable_waddons[,plotmetrics])
    # Add in NotAvailable instead of NA to make sorting better
    labeltable <- do.call(cbind.data.frame, lapply(labeltable_temp, function(x) {
        factor(ifelse(is.na(x), "NotAvailable", as.character(x)), levels = c("NotAvailable", levels(x)))
    }))
    rownames(labeltable) <- rownames(labeltable_temp)
    labeltable <- labeltable[rev(do.call(order, labeltable)),]
    
    customcolorlist <- custom_annotation_list_from_colorguide(COI = plotmetrics, colorguide)
    annotationlist1 <- annotationlist_builder(labeltable, customcolorlist = customcolorlist)
    ## Dummy table to work with function.
    dummytab <- labeltable
    dummytab[] <- 1
    labelplot <- create_heatmap(counttab = dummytab, subsetnum = FALSE, scale_data = FALSE, colmetatable = NULL, colannotationlist = NULL, 
                                rowmetatable = labeltable, rowannotationlist = annotationlist1, 
                                colclusterparam = FALSE, rowclusterparam = FALSE, separate_legend = FALSE, heatmapcolorparam = NULL)
    return(labelplot)
    
}


# --------------------------------- In group vs Out group enrichment tests --------------------------------- 

ingroup_vs_outgroup_cohort_enrichment_tests <- function(tabulation_table, group_column, cohort_tabulations_outpath) {
    summary_table_pval_outlist <- summary_table_sign_outlist <- summary_table_ratio_outlist <- summary_ci_table_outlist <- list()
    for (nmfnum in seq_len(length(sort(unique(na.omit(tabulation_table[,group_column])))))) {
        nmfsel <- sort(unique(na.omit(tabulation_table[,group_column])))[nmfnum]
        # tabulation_table_sel <- tabulation_table[SOIrna,]
        tabulation_table_sel <- tabulation_table
        tabulation_table_sel[,group_column] <- ifelse(tabulation_table_sel[,group_column] == nmfsel, nmfsel, paste0("NOT_", nmfsel))
        
        summarystatfile <- paste0(cohort_tabulations_outpath, "COI_tabulation_", "group_column", "_", nmfsel, ".csv")
        sumstattab_seq <- summarize_table(intable = tabulation_table_sel, groupvar = group_column, 
                                          outfile = summarystatfile, calc_stats = TRUE, calc_ci = TRUE)
        
        collapsed_summary_table <- cbind(category = rep(names(sumstattab_seq), lapply(sumstattab_seq, nrow)), do.call(rbind.fill, sumstattab_seq))
        collapsed_summary_table[collapsed_summary_table[,paste0(group_column, "_cat")] == "Min.", paste0(group_column, "_cat")] <- "NUMERIC"
        collapsed_summary_table[,"combined_catval"] <- apply(collapsed_summary_table[,c("category", paste0(group_column, "_cat"))], 1, 
                                                             function(x) paste(x, collapse = "__"))
        ## First one - is a pval table
        summary_table_pval <- collapsed_summary_table[!collapsed_summary_table[,"statval"] %in% c(NA, ""),c("category", paste0(group_column, "_cat"), "statval", "combined_catval")]
        summary_table_pval_outlist[[nmfnum]] <- cbind(summary_table_pval[,c("combined_catval", "statval")], cluster = nmfsel)
        names(summary_table_pval_outlist)[nmfnum] <- nmfsel
        
        ## Second one - is a CI table
        summary_table_ci_val <- collapsed_summary_table[!collapsed_summary_table[,"ci_val"] %in% c(NA, ""),c("category", paste0(group_column, "_cat"), "ci_val", "combined_catval")]
        summary_ci_table_outlist[[nmfnum]] <- cbind(summary_table_ci_val[,c("combined_catval", "ci_val")], cluster = nmfsel)
        names(summary_ci_table_outlist)[nmfnum] <- nmfsel
        
        ## For each category - 
        signoutlist <- ratiooutlist <- list()
        for (combined_catval_num in seq_len(length(summary_table_pval[,"combined_catval"]))) {
            combined_catval_sel <- summary_table_pval[,"combined_catval"][combined_catval_num]
            if (grepl("__NUMERIC", combined_catval_sel)) {
                catsel <- gsub("NUMERIC", "Mean", combined_catval_sel)
                datasel <- collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == catsel,]
                signout <- ifelse(datasel[,nmfsel] > datasel[,paste0("NOT_", nmfsel)], 1, -1)
                
                # ratioout <- log2(as.numeric(datasel[,nmfsel]) / as.numeric(datasel[,paste0("NOT_", nmfsel)]))
                # if (is.infinite(ratioout)) {ratioout <- 10 * sign(ratioout)}
                ratioout <- as.numeric(datasel[,nmfsel]) / as.numeric(datasel[,paste0("NOT_", nmfsel)])
                if (is.infinite(ratioout)) {ratioout <- 10}
                
            } else {
                datasel <- collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == combined_catval_sel,]
                outratios <- as.numeric(unlist(datasel[,c(nmfsel, paste0("NOT_", nmfsel))])) / 
                    # as.numeric(unlist(collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == "PATNUM__total",
                    #                                           c(nmfsel, paste0("NOT_", nmfsel))]))
                    as.numeric(unlist(collapsed_summary_table[grepl("__total", collapsed_summary_table[,"combined_catval"]),
                                                              c(nmfsel, paste0("NOT_", nmfsel))]))
                signout <- ifelse(outratios[1] > outratios[2], 1, -1)
                
                # ratioout <- log2(outratios[1] / outratios[2])
                # if (is.infinite(ratioout)) {ratioout <- 10 * sign(ratioout)}
                ratioout <- outratios[1] / outratios[2]
                if (is.infinite(ratioout)) {ratioout <- 10}
                
            }
            signoutlist[[combined_catval_num]] <- data.frame(combined_catval = combined_catval_sel, cluster = nmfsel, statsign = as.numeric(signout))
            ratiooutlist[[combined_catval_num]] <- data.frame(combined_catval = combined_catval_sel, cluster = nmfsel, ratio = as.numeric(ratioout))
            names(signoutlist)[combined_catval_num] <- names(ratiooutlist)[combined_catval_num] <- combined_catval_sel
        }
        signouttable <- do.call(rbind, signoutlist)
        summary_table_sign_outlist[[nmfnum]] <- signouttable
        names(summary_table_sign_outlist)[nmfnum] <- nmfsel
        
        ratioouttable <- do.call(rbind, ratiooutlist)
        summary_table_ratio_outlist[[nmfnum]] <- ratioouttable
        names(summary_table_ratio_outlist)[nmfnum] <- nmfsel
        
        # summary_table_value <- collapsed_summary_table[!summary_table_pval[,"statval"] %in% c(NA, ""),c("category", "rnacluster_cat", "statval")]
        # summary_table_pval[summary_table_pval[,"rnacluster_cat"] == "Min.","rnacluster_cat"] <- "NUMERIC"
    }
    # Summary sign out table
    summary_table_sign_temp <- do.call(rbind, summary_table_sign_outlist)
    summary_table_sign_temp2 <- dcast(summary_table_sign_temp, combined_catval ~ cluster, value.var = "statsign")
    summary_table_sign_FULL <- apply(summary_table_sign_temp2[,2:ncol(summary_table_sign_temp2)], 2, function(x) as.numeric(as.character(x)))
    rownames(summary_table_sign_FULL) <- summary_table_sign_temp2[,"combined_catval"]
    summary_table_sign_FULL <- summary_table_sign_FULL[summary_table_pval_outlist[[1]][,1], ]
    
    # ratios instead of signs
    summary_table_ratio_temp <- do.call(rbind, summary_table_ratio_outlist)
    summary_table_ratio_temp2 <- dcast(summary_table_ratio_temp, combined_catval ~ cluster, value.var = "ratio")
    summary_table_ratio_FULL <- apply(summary_table_ratio_temp2[,2:ncol(summary_table_ratio_temp2)], 2, function(x) as.numeric(as.character(x)))
    rownames(summary_table_ratio_FULL) <- summary_table_ratio_temp2[,"combined_catval"]
    summary_table_ratio_FULL <- summary_table_ratio_FULL[summary_table_pval_outlist[[1]][,1], ]
    
    summary_table_pval_temp <- do.call(rbind, summary_table_pval_outlist)
    summary_table_pval_temp2 <- dcast(summary_table_pval_temp, combined_catval ~ cluster, value.var = "statval")
    summary_table_pval_FULL <- apply(summary_table_pval_temp2[,2:ncol(summary_table_pval_temp2)], 2, function(x) as.numeric(as.character(x)))
    rownames(summary_table_pval_FULL) <- summary_table_pval_temp2[,"combined_catval"]
    summary_table_pval_FULL <- summary_table_pval_FULL[summary_table_pval_outlist[[1]][,1], ]
    
    summary_table_cival_temp <- do.call(rbind, summary_ci_table_outlist)
    summary_table_cival_temp2 <- dcast(summary_table_cival_temp, combined_catval ~ cluster, value.var = "ci_val")
    summary_table_cival_FULL <- apply(summary_table_cival_temp2[,2:ncol(summary_table_cival_temp2)], 2, function(x) as.character(x))
    rownames(summary_table_cival_FULL) <- summary_table_cival_temp2[,"combined_catval"]
    summary_table_cival_FULL <- summary_table_cival_FULL[summary_ci_table_outlist[[1]][,1], ]
    
    # Return just the CI values (no mean)
    summary_table_cival_FULL_out <- summary_table_cival_FULL
    summary_table_cival_FULL_out[] <- apply(summary_table_cival_FULL_out, c(1,2), function(x) sub(".*? ", "", x))
    
    # summary_table_pval_FULL_signed <- summary_table_pval_FULL * summary_table_sign_FULL
    summary_table_pval_FULL_signed <- summary_table_ratio_FULL
    
    # summary_table_pval_FULL_filtered <- summary_table_pval_FULL
    summary_table_pval_FULL_filtered_signed <- summary_table_ratio_FULL
    summary_table_pval_FULL_filtered_signed[summary_table_pval_FULL > 0.05] <- NA
    # summary_table_pval_FULL_filtered_signed <- summary_table_pval_FULL_filtered * summary_table_sign_FULL
    
    return(list(summary_table_pval_FULL_filtered_signed = summary_table_pval_FULL_filtered_signed,
                summary_table_cival_FULL_out = summary_table_pval_FULL_filtered_signed,
                summary_table_pval_FULL_signed = summary_table_pval_FULL_signed))
    
}

# --------------------------------- hazard ratio heatmap --------------------------------- 

hr_heatmap <- function(hr_table) {
    ## Lets output a cleaned version of the cumhazouttable, and a heatmap with the CI and text within it
    # tt1 <- data.frame(fullcumhazouttable)
    genestoagg <- unique(gsub("_.*", "", rownames(fullcumhazouttable)))
    ## I can definitely aggregate this, but idk how, so I am gonna use a forloop instead.
    HRoutlist <- textoutlist <- list()
    for (genenum in seq_len(length(genestoagg))) {
        ## Grab the rows we need
        genesel <- genestoagg[genenum]
        selrows <- fullcumhazouttable[grepl(genesel, rownames(fullcumhazouttable)),]
        
        ## Output to a new table
        # textoutlist[[genenum]] <- apply(selrows, 2, function(x){paste0(round(x[1], 2), " (",round(x[2], 2), " - ", round(x[3], 2), ")")})
        ## Lets try just the CIs
        textoutlist[[genenum]] <- apply(selrows, 2, function(x){paste0("(",round(x[2], 2), " - ", round(x[3], 2), ")")})
        names(textoutlist)[genenum] <- genesel
        
        ## Write out just the HR to a list
        HRoutlist[[genenum]] <- selrows[grepl("_exp", rownames(selrows)),]
        names(HRoutlist)[genenum] <- genesel
    }
    HRtable <- hr_table
    HRtexttable <- do.call(rbind, textoutlist)
    
    heatmapcolorparam <- colorRamp2(c(0, 1, 6), c( "#0000A7", "White", "#A30000"))
    hrmap = Heatmap(HRtable[rownames(maptab),], 
                    border = TRUE,
                    rect_gp = gpar(col = "black", lwd = 1),
                    col = heatmapcolorparam,    ## Define the color scale for the heatmap
                    row_title = "Gene",                                       ## Name the rows
                    column_title = "Event",                                  ## Name the columns
                    
                    cluster_columns = FALSE,                         ## Cluster the columns or leave as is
                    cluster_rows = FALSE,                            ## Cluster the rows or leave as is
                    
                    show_column_names = TRUE  ,                               ## Show the Column Names
                    column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
                    show_row_names = TRUE,                                   ## Show the row names
                    row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
                    row_names_gp = gpar(fontsize=6),
                    
                    heatmap_legend_param = list(
                        title = "Hazard Ratio",
                        legend_height = unit(2.5, "cm"),
                        title_gp = gpar(fontsize = 8, fontface = "bold")),
                    
                    height = unit(20,"cm"),
                    width = unit(10,"cm"),
                    
                    ## Adds in the cor value from the maptab
                    # cell_fun = function(j, i, x, y, width, height, fill) {
                    #     grid.text(paste0(sprintf("%.2f", t(moduleTraitCor)[,colnames(maptab),drop=FALSE][i, j]), "\n", "(",
                    #                      sprintf("%.2f", t(moduleTraitPvalue)[,colnames(maptab),drop=FALSE][i, j]), ")"),
                    #               x, y, gp = gpar(fontsize = 8))
                    # }
                    cell_fun = function(j, i, x, y, width, height, fill) {
                        # grid.text(paste0(sprintf("%.2f", moduleTraitCor[i, j]), "\n", "(",
                        #                  sprintf("%.2f", moduleTraitPvalue[i, j]), ")"),
                        grid.text(HRtexttable[rownames(maptab),][i, j],
                                  x, y, gp = gpar(fontsize = 6))
                    }
                    
    )
    hmoutfile2 <- paste0(outfilepathsurvival, "HR_heatmap.pdf")
    # hmoutfile <- paste0(outfilepathsurvival, "CORRECTED_event_heatmap_log2fcsort_with_annot.pdf")
    pdf(hmoutfile2, 11, 9)
    draw(hrmap)
    junk <- dev.off()
    
}



# --------------------------------- cluster module meta summary heatmap --------------------------------- 

cluster_module_meta_summary_heatmap_function <- function(ingenestocolortable, incounttable, bioreptable_waddons, 
                                                        selected_clusterlabel, eigengeneorder, metacols_sel, ...) {
    ## Read in eigengene table
    # GOI <- rownames(ingenestocolortable[ingenestocolortable[,"moduleColors"] != "grey",]) ## Only plotting on genes that are NOT grey
    GOI <- rownames(ingenestocolortable[ingenestocolortable[,"moduleColors"] %in% eigengeneorder,]) ## Only plotting on genes that are NOT grey
    GOI <- GOI[GOI %in% rownames(incounttable)]
    plotSOI <- colnames(incounttable)
    
    # Now lets make a summary - where we have a split heatmap - with modules labeled and split on rows, and then sample labeled and split on columns
    hmplottab <- incounttable[GOI, plotSOI]
    colsplitparam <- factor(bioreptable_waddons[plotSOI, selected_clusterlabel])
    rowsplitparam <- factor(ingenestocolortable[GOI,"moduleColors"], levels = eigengeneorder)
    
    sampleannottable <- merge(bioreptable_waddons[plotSOI, metacols_sel], 
                              bioreptable_waddons[plotSOI, selected_clusterlabel, drop=FALSE], by = "row.names")
    rownames(sampleannottable) <- sampleannottable[,"Row.names"]
    sampleannottable <- sampleannottable[plotSOI,!grepl("Row.names", colnames(sampleannottable))]
    
    ## Make the annotation (the main plot we want)
    sampleannotationlist <- custom_annotation_list_from_colorguide(COI = colnames(sampleannottable), colorguide)
    
    ## Top annotation
    temp1 <- vector("list", length(sampleannotationlist))
    names(temp1) = names(sampleannotationlist)
    annotlegendlist = lapply(temp1, function(x) x[[1]] =
                                 list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
    
    # ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
    for (annotnum in seq_len(length(annotlegendlist))) {
        annotname_select <- names(annotlegendlist)[annotnum]
        if (!is.null(names(sampleannotationlist[[annotname_select]]))) {
            annotlegendlist[[annotnum]] <- c(annotlegendlist[[annotnum]], list(at = names(sampleannotationlist[[annotname_select]])))
        } else {next}
    }
    ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
    
    ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
    showlegendparam = unname(unlist(lapply(sampleannotationlist, function(x) {
        numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
        is.null(numterms) || numterms <= 10})))
    hatop = HeatmapAnnotation(df = sampleannottable,
                              col = sampleannotationlist,
                              na_col = "white",
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                              annotation_name_side = "left",
                              simple_anno_size = unit(min(60/length(sampleannotationlist), 5),"mm"),
                              show_legend = TRUE,
                              annotation_legend_param = annotlegendlist)
    
    # gene annotation
    geneannottable <- ingenestocolortable[GOI,"moduleColors", drop=FALSE]
    genecustomcolorlist <- list(moduleColors = eigengeneorder)
    names(genecustomcolorlist[["moduleColors"]]) <- eigengeneorder
    geneannotationlist <- annotationlist_builder(geneannottable, customcolorlist = genecustomcolorlist)
    
    ## Side annotation
    ## Define parameters for each of the labels on the annotation bars
    temp1 <- vector("list", length(geneannotationlist))
    names(temp1) = names(geneannotationlist)
    annotlegendlist = lapply(temp1, function(x)
        x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))
    
    # ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
    for (annotnum in seq_len(length(annotlegendlist))) {
        annotname_select <- names(annotlegendlist)[annotnum]
        if (!is.null(names(geneannotationlist[[annotname_select]]))) {
            annotlegendlist[[annotnum]] <- c(annotlegendlist[[annotnum]], list(at = names(geneannotationlist[[annotname_select]])))
        } else {next}
    }
    ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
    
    ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
    showlegendparam = unname(unlist(lapply(geneannotationlist, function(x) {
        numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
        is.null(numterms) || numterms <= 10})))
    
    ## Look for any empty annotations - fill them with white, and later, make sure to hide their legend
    emptyannots = names(sapply(geneannotationlist, length)[sapply(geneannotationlist, length)==0])
    if (length(emptyannots) > 0){
        for (i in 1:length(emptyannots)) {
            temp1 = "white"
            names(temp1) = emptyannots[i]
            geneannotationlist[[emptyannots[i]]] = temp1
        }
        showlegendparam[which(names(geneannotationlist) %in% emptyannots)] = FALSE
    }
    ## Add param that will bolden the side annotation bars if it is <100, and omit the grid lines if more
    if (nrow(geneannottable) < 100) {
        sideannotation_linebold_param <- gpar(fontsize = 0.5)} else {sideannotation_linebold_param <- NULL}
    haside = rowAnnotation(df = geneannottable,
                           col = geneannotationlist,
                           na_col = "white",
                           # gp = gpar(fontsize = 0.01),
                           gp = sideannotation_linebold_param,
                           show_annotation_name=TRUE,
                           annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                           annotation_name_side = "top",
                           simple_anno_size = unit(min(60/length(geneannotationlist), 5),"mm"),
                           show_legend = TRUE,
                           annotation_legend_param = annotlegendlist)
    
    
    hmplottab = as.matrix(t(apply(hmplottab, 1, function(x) zscore(x))))
    heatmapcolorparam <- colorRamp2(breaks = c(4, 0, -4), c("darkred", "white", "darkblue"))
    summary_outhm <- Heatmap(as.matrix(hmplottab),
                             col = heatmapcolorparam,
                             row_title = "Genes",column_title = "Samples",
                             row_split = rowsplitparam, column_split = colsplitparam,
                             top_annotation = hatop,
                             
                             cluster_columns = TRUE, cluster_rows = TRUE,
                             cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                             
                             show_row_dend = FALSE, show_column_dend = FALSE,
                             
                             # rect_gp = gpar(col = "black", lwd = 0.5),
                             
                             show_column_names = FALSE, column_names_gp = gpar(fontsize = 6), 
                             show_row_names = FALSE,
                             
                             heatmap_legend_param = list(
                                 title = "Zscore",
                                 title_gp = gpar(fontsize = 8, fontface = "bold")),
                             height = unit(min((nrow(hmplottab)/2), 12),"cm"),
                             width = unit(min(ncol(hmplottab), 18),"cm")
                             
    )
    
    out1 <- draw(summary_outhm + haside, annotation_legend_side = "bottom", padding = unit(c(5, 20, 5, 5), "mm"))
    return(out1)
    # pdf(paste0(outfilepathintegration, "nmf_eigen_meta_summary_heatmap.pdf"), useDingbats = FALSE, width = 20, height = 15)
    # draw(summary_outhm + haside)
    # junk <- dev.off()
}


# --------------------------------- eigengene vs cluster heatmap --------------------------------- 

eigen_v_cluster_heatmap <- function(eigengenecounttable, selected_clusterlabel, eigengeneorder, colorguide) {
    
    # this is the plottable that we need
    hmplottable <- t(eigengenecounttable[,paste0("ME", eigengeneorder)])
    plotSOI <- colnames(hmplottable)
    
    # But i also want the rank table by cluster, so i need a box plot style melted table, and then go from there
    temptab <- na.omit(merge(bioreptable_waddons[,selected_clusterlabel, drop = FALSE], eigengenecounttable, by = "row.names"))
    rownames(temptab) <- temptab[,"Row.names"]
    temptab <- temptab[,!grepl("Row.names", colnames(temptab))]
    cluster_eigenrank_table_value <- aggregate(temptab[,colnames(eigengenecounttable)], by = list(temptab[,selected_clusterlabel]), mean)
    rownames(cluster_eigenrank_table_value) <- cluster_eigenrank_table_value[,"Group.1"]
    cluster_eigenrank_table_value <- cluster_eigenrank_table_value[,!grepl("Group.1", colnames(cluster_eigenrank_table_value))]
    
    # I think I want to output a helper table here - that will give relative values between them
    ## Namely, I want these as scaled (???) do I?? values, but most importantly, I want the difference between the highest value, and the NEXT highest value ...
    max_minus_next_diff_table <- data.frame(t(apply(cluster_eigenrank_table_value, 2, function(x) {
        # x <- cluster_eigenrank_table_value[,1]
        # rankvals <- rank(x)
        max_minus_next_diff <- x[which.max(x)] -  x[-which.max(x)][which.max(x[-which.max(x)])]
        maxval <- paste0("meth_3cluster_", which.max(x))
        c(maxval, max_minus_next_diff)
    })))
    colnames(max_minus_next_diff_table) <- c("maxcluster", "max_minus_next_diff")
    max_minus_next_diff_table <- max_minus_next_diff_table[paste0("ME", eigengeneorder),,drop=FALSE]
    max_minus_next_diff_table <- max_minus_next_diff_table[rev(order(max_minus_next_diff_table[,"maxcluster"], 
                                                                 -as.numeric(max_minus_next_diff_table[,"max_minus_next_diff"]), decreasing = TRUE)),,drop=FALSE]
    # max_minus_next_diff_table
    # max_minus_next_diff_table <- max_minus_next_diff_table[order(max_minus_next_diff_table[,"max_minus_next_diff"], decreasing = TRUE),,drop=FALSE]
    # max_minus_next_diff_table <- max_minus_next_diff_table[order(max_minus_next_diff_table[,"max_minus_next_diff"], decreasing = TRUE),]
    
    
    cluster_eigenrank_table_rank <- data.frame(apply(cluster_eigenrank_table_value, 2, function(x) rev(rank(x))))
    # Gives an attempt at ordering - will probably need to be manually ordered though
    cluster_eigenrank_table_rank <- data.frame(t(cluster_eigenrank_table_rank[,do.call(order, data.frame(t(cluster_eigenrank_table_rank)))]))
    
    # col annotation is just the cluster labels (easy)
    colmetatable <- bioreptable_waddons[plotSOI, selected_clusterlabel, drop = FALSE]
    colsplitparam <- data.frame(factor(bioreptable_waddons[plotSOI, selected_clusterlabel]))
    dimnames(colsplitparam) <- list(plotSOI, selected_clusterlabel)
    colannotationlist <- custom_annotation_list_from_colorguide(COI = selected_clusterlabel, colorguide)
    
    # USING THIS RANK TABLE for annotation
    # Black to white color range
    rowmetatable <- cluster_eigenrank_table_rank[paste0("ME", eigengeneorder), names(colannotationlist[[selected_clusterlabel]])]
    blackwhite_color_range <- colorRampPalette(c("#e5e5e5", "#191919"))
    rowcustomcolorlist <- rep(list(colorRamp2(breaks = sort(unique(unlist(rowmetatable))),
                                              # colors = brewer.pal(length(unique(unlist(rowmetatable))), "Greys")
                                              colors = blackwhite_color_range(length(unique(unlist(rowmetatable)))) 
                                              )),
                              length(unique(unlist(rowmetatable))))
    names(rowcustomcolorlist) <- colnames(rowmetatable)
    rowannotationlist <- annotationlist_builder(rowmetatable, customcolorlist = rowcustomcolorlist)
    
    #' customcolorlist = list(A = c("high" = "red", "low" = "blue"),
    #'                        #'                        B = circlize::colorRamp2(seq(-5, 5, length = 3),
    #'                        #'                        RColorBrewer::brewer.pal(3, "Reds")))
    
    heatmapcolorparam <- colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
    outhm1 <- create_heatmap(counttab = hmplottable, scale_data = FALSE, separate_legend = FALSE,
                             colmetatable = colmetatable, colannotationlist = colannotationlist,
                             rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                             colclusterparam = TRUE, rowclusterparam = FALSE,
                             heatmapcolorparam = heatmapcolorparam,
                             columnsplitparam = colsplitparam)
    return(list(heatmap = outhm1[[1]], cluster_eigenrank_table_rank = rowmetatable,
                max_minus_next_diff_table = max_minus_next_diff_table, cluster_eigenrank_table_value = cluster_eigenrank_table_value))
    
}



# --------------------------------- create boxplots for every feature by annotation column --------------------------------- 

# Eigengene boxplots across cluster ... really boxplots for any annotation with a count table
## So input is a counttable, and each feature will get a boxplot, with the separate boxplots being pulled from an annotation table by column
## output will be N boxplot figures (N = features) by X folders (X = number of annotation columns)
## Also output the spearman correlation plot and values as well
## Need the outfilepath - cause we write out a bajillion boxplots, and dont want to save those all to and output object

boxplots_from_counttable_by_annotation <- function(counttable, boxplot_annotationtable, outfilepath, 
                                                   calculate_corrvalues = FALSE, 
                                                   bp_color_ref = NULL) {
    # These are the features and annotations we will run over
    bp_features <- rownames(counttable)
    bp_annots <- colnames(boxplot_annotationtable)
    
    params_grid <- expand.grid(bp_features, bp_annots)
    colnames(params_grid) <- c("bp_features", "bp_annots")
    wilcox_statoutlist <- spearman_statoutlist <- list()
    for (params_num in seq_len(nrow(params_grid))) {
    # for (params_num in 1:75) {
        
        # Grab this runs parameters
        param_sel <- params_grid[params_num,]
        feature_sel <- param_sel[,"bp_features"]
        annot_sel <- param_sel[,"bp_annots"]
        dir.create(paste0(outfilepath, annot_sel, "/"), showWarnings = FALSE, recursive = TRUE)
        
        counttable_sel <- t(counttable[feature_sel,,drop=FALSE])
        annottable_sel <- boxplot_annotationtable[,annot_sel,drop=FALSE]
        # nmf_sel <- nmf_clustermembership_table[,"combined",drop=FALSE]
        bptab <- na.omit(merge(counttable_sel, annottable_sel, by = "row.names")[,c(3,1,2)])
        bptab[,2] <- c("cohort")

        ## IF THE FIRST COL IS CONT (like for labs) - then split into quartiles
        rawspearmanstatout <- NULL
        if (is.numeric(bptab[,1])) {
            ## Also sneak in an actual corr with the raw values here jsut to see
            scattertab <- bptab[,c(1,3)]
            scattertab[,3] <- gsub("ME", "", colnames(scattertab)[2])

            scatterout <- scatter_plotter(indata = scattertab[,c(1,2)], colorvar = scattertab[,3,drop=FALSE],
                                          labsparam = list(title = paste0(feature_sel, " values for samples by ", annot_sel), x = annot_sel, y = feature_sel),
                                          plotstats = TRUE)
            pdf(paste0(outfilepath, annot_sel, "/", feature_sel, "_rawval_scatter.pdf"), useDingbats = FALSE)
            print(scatterout)
            junk <- dev.off()

            ## And capture this spearman value
            rawspearmanstatval <- cor.test(scattertab[,1], scattertab[,2], method = "spearman", exact = FALSE)
            rawspearmanstatout <- c(paste0(annot_sel, "_raw"), feature_sel, rawspearmanstatval$p.value, rawspearmanstatval$estimate)

            ## Create the bptab
            cattab <- cut2(t(bptab[,1]),g=4)
            levels(cattab)[match(levels(cattab)[1],levels(cattab))] <- paste0(annot_sel, "_Q1min")
            levels(cattab)[match(levels(cattab)[2],levels(cattab))] <- paste0(annot_sel, "_Q2")
            levels(cattab)[match(levels(cattab)[3],levels(cattab))] <- paste0(annot_sel, "_Q3")
            levels(cattab)[match(levels(cattab)[4],levels(cattab))] <- paste0(annot_sel, "_Q4max")
            bptab[,1] <- as.character(cattab)
        }

        ## Pre calc the avg per group so we can sort in that order:
        avgfeatvalue <- aggregate(bptab[,3], by = list(bptab[,1]), mean)
        avgfeatvalue <- avgfeatvalue[order(avgfeatvalue[,2]),]
        avgfeatvalue[,"Eigen_Score_Rank"] <- seq(1:nrow(avgfeatvalue))
        colnames(avgfeatvalue) <- c("Group", "featurevalue", "group_feature_rank")
        
        # Add custom coloring for this particular section
        # bp_color_ref <- custom_annotation_list_from_colorguide(COI=group_column_selected, colorguide)
        if (!is.null(bp_color_ref)){
            if (annot_sel %in% names(bp_color_ref)) {
                featorder_param <- names(bp_color_ref[[annot_sel]])
            }
        } else {
            featorder_param <- levels(bptab[,1])
        }
        # Add custom coloring for this particular section

        ## OUTPUT THIS
        bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "category",
                                 labsparam = list(title = paste0(feature_sel, " values for samples by ", annot_sel), x = annot_sel, y = feature_sel,
                                                  catorder = "cohort", featorder = featorder_param),
                                 plotstats = "intra", testtypeparam = "wilcox.test"
        )
        
        # Add custom coloring for this particular section
        # bp_color_ref <- custom_annotation_list_from_colorguide(COI=group_column_selected, colorguide)
        if (!is.null(bp_color_ref)){
            if (annot_sel %in% names(bp_color_ref)) {
                bpout <- bpout + scale_fill_manual(breaks = names(bp_color_ref[[annot_sel]]), values = unname(bp_color_ref[[annot_sel]]))
            }
        }
        # Add custom coloring for this particular section
        
        pdf(paste0(outfilepath, annot_sel, "/", feature_sel, "_boxplot.pdf"), useDingbats = FALSE)
        print(bpout)
        junk <- dev.off()
        # ## Save out the bptab
        # bptablist[[params_num]] <- bptab

        ## Now we also want a corr plot if applicable?
        if (calculate_corrvalues) {
            scattertab <- bptab[,c(1,3)]
            scattertab[,3] <- avgfeatvalue[,3][match(scattertab[,1], avgfeatvalue[,1])]
            scattertab[,4] <- colnames(scattertab)[2]
            
            scatterlabelparam <- as.vector(avgfeatvalue[,1])
            names(scatterlabelparam) <- avgfeatvalue[,3]
            scatterout <- scatter_plotter(indata = scattertab[,c(3,2)], colorvar = scattertab[,4,drop=FALSE],
                                          labsparam = list(title = paste0(feature_sel, " values for samples by ", annot_sel), x = annot_sel, y = feature_sel),
                                          plotstats = TRUE, addjitterparam = 0.1)
            scatterout <- scatterout + scale_x_continuous(breaks = avgfeatvalue[,3], labels = avgfeatvalue[,1])
            # scatterout <- scatterout + geom_point(alpha = 0.7, shape = 16, position = position_jitter(width = 0.1))
            pdf(paste0(outfilepath, annot_sel, "/", feature_sel, "_scatter.pdf"), useDingbats = FALSE)
            print(scatterout)
            junk <- dev.off()
            
            # Repeat for the spearman correlation
            spearmanstatout <- cor.test(scattertab[,3], scattertab[,2], method = "spearman", exact = FALSE)
            if (!is.null(rawspearmanstatout)) { ## Attach the raw spearman corr out as well if applicable
                spearman_statoutlist[[params_num]] <- rbind(rawspearmanstatout, c(paste0(annot_sel, "_quartile"), feature_sel, spearmanstatout$p.value, spearmanstatout$estimate))
            } else {
                spearman_statoutlist[[params_num]] <- c(as.character(annot_sel), as.character(feature_sel),
                                                        spearmanstatout$p.value, spearmanstatout$estimate)
            }
            names(spearman_statoutlist)[params_num] <- feature_sel
        }
        

        # Grab the stats into a table
        # Create the combinatorial table for all combos of the groups
        # If this is a factor - then keep the levels of the factor for the order 2022-09-28
        if (is.factor(bptab[,1])) {
            combtab = combn(na.omit(levels(bptab[,1])), 2, simplify=F)
        } else {
            combtab = combn(as.character(unique(bptab[!is.na(bptab[,1]),1])), 2, simplify=F)
        }
        
        
        # Apply a functiona cross each combo using the bptab and pulling out each group
        wilcoxstatsout <- lapply(combtab, function(x){
            group1 <- bptab[bptab[,1] %in% x[1], 3]
            group2 <- bptab[bptab[,1] %in% x[2], 3]
            out1 <- c(as.character(annot_sel), x[1], x[2], as.character(feature_sel), wilcox.test(group1, group2)$p.value)
            out2 <- cohens_d(group1, group2)
            out3 <- c(mean(group1, na.rm = TRUE), mean(group2, na.rm = TRUE))
            c(out1, out2, out3)
        })
        wilcox_statoutlist[[params_num]] <- do.call(rbind, wilcoxstatsout)

        ## Name all of our outlists
        names(wilcox_statoutlist)[params_num] <- feature_sel

    }
    # # Cat our lists
    wilcox_summarytable <- do.call(rbind, wilcox_statoutlist)
    colnames(wilcox_summarytable) <- c("Cohort", "Feature1", "Feature2", "Feature", "wilcox_pval", "cohens_d", "Feature1_mean", "Feature2_mean")
    if (calculate_corrvalues) {
        spearman_summarytable <- do.call(rbind, spearman_statoutlist)
        colnames(spearman_summarytable) <- c("Cohort", "Feature", "spearman_pval", "spearman_rval")
    } else {spearman_summarytable <- NULL}
    
    return(list(wilcox_summarytable = wilcox_summarytable, spearman_summarytable = spearman_summarytable))
    
    # eigenrank_summarytable <- do.call(rbind, eigenranklist)
    # eigenrank_summarytable[,"eigengene"] <- gsub("\\..*", "", rownames(eigenrank_summarytable))
    # eigenrank_summarytable <- eigenrank_summarytable[order(eigenrank_summarytable[,"Group"], eigenrank_summarytable[,"Eigen_Score_Rank"], decreasing = TRUE),]
    # 
    # # Write out the results
    # write.table(wilcox_summarytable, paste0(COIbpoutpath, COItab_label,"_wilcox_summary_table.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    # write.table(spearman_summarytable, paste0(COIbpoutpath, COItab_label, "_spearman_summary_table.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    # write.table(eigenrank_summarytable, paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_table.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    # 
    # #     bptablist[[params_num]] <- bptab
    # 
    # ALL_wilcox_statoutlist[[COItabnum]] <- wilcox_summarytable
    # names(ALL_wilcox_statoutlist[COItabnum]) <- COItab_label
    # ALL_spearman_statoutlist[[COItabnum]] <- spearman_summarytable
    # names(ALL_spearman_statoutlist[COItabnum]) <- COItab_label
    # 
    # }
    # ALL_wilcox_sumarytable <- do.call(rbind, ALL_wilcox_statoutlist)
    # ALL_spearman_sumarytable <- do.call(rbind, ALL_spearman_statoutlist)
    
}


# --------------------------------- cluster vs event heatmap --------------------------------- 


# what do i want to control - selected_cluster_table, events (EOI), cohorts (COI), outfilepath, bioreptable (indexed with EOI and COI)
# but i need new E and C ADDED TO the bioreptable already

cluster_vs_event_kmanalysis_summary_heatmap <- function(selected_cluster_table, COIref_tablist, survivalplot_outfilepath, bioreptable_waddons,
                                                        eventpvalcutoff, return_output_tables = FALSE, EOI = NULL, ...) {

    ## Grab EOI from our bioreptable
    if (!is.null(EOI)) {EOI <- gsub("^C_", "", colnames(bioreptable_waddons)[grepl("^C_", colnames(bioreptable_waddons))])}
    EOIcols <- apply(expand.grid(c("T_", "C_"), EOI, stringsAsFactors = FALSE), 1, paste, collapse = "")

    # For each cluster label - lets run this analysis
    expanded_nmflabeltable <- data.frame(make_comparison_columns(selected_cluster_table[,2,drop=FALSE]))
    # Turn our expanded nmfcluster table back into a list with each column as an entry:
    expanded_nmflabeltable_list <- lapply(expanded_nmflabeltable[,!colnames(expanded_nmflabeltable) %in% "PATNUM"], function(x)
        data.frame(x, row.names = rownames(expanded_nmflabeltable[,!colnames(expanded_nmflabeltable) %in% "PATNUM"])))
    for (clusterlabelnum in seq_len(length(colnames(expanded_nmflabeltable[,!colnames(expanded_nmflabeltable) %in% "PATNUM"])))){
        colnames(expanded_nmflabeltable_list[[clusterlabelnum]]) <- 
            colnames(expanded_nmflabeltable[,!colnames(expanded_nmflabeltable) %in% "PATNUM"])[clusterlabelnum]
    }

    ## Ok - so we want to iterate over every GOI, and every event, and see what comes out
    survival_fulltable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE), list(
        bioreptable_waddons, cbind.data.frame(PATNUM = rownames(expanded_nmflabeltable), expanded_nmflabeltable)
    ))
    write.table(survival_fulltable, paste0(survivalplot_outfilepath, "full_survival_intable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    ## COI tab list - the categories to compare against
    COItablist <- c(COIref_tablist, expanded_nmflabeltable_list)

    fullpvaloutlist <- fullcumhazoutlist <- fullcohortNtablelist <- list()
    clusterplotCOI <- c()
    for (COInum in seq_len(length(COItablist))) {
        ## Select COI
        COIsel <- COItablist[[COInum]]
        COIlabel <- names(COItablist)[COInum]

        # Create the intable with events and the specific category of interest
        survival_intable <- merge(survival_fulltable[,c("PATNUM", EOIcols)], cbind(PATNUM = rownames(COIsel), COIsel), by = "PATNUM")
        colnames(survival_intable)[ncol(survival_intable)] <- COIlabel

        survpvallist <- survHRlist <- cohortNtablelist <- list()
        for (eventnum in seq_len(length(EOI))) {
            ## Select event
            eventsel <- EOI[eventnum]

            ## Create the subsurvivaltab
            survivaldata <- na.omit(data.frame(survival_intable[,c(paste0(c("T_", "C_"), eventsel), COIlabel, "PATNUM")]))
            rownames(survivaldata) <- survivaldata[,"PATNUM"]
            survivaldata <- survivaldata[,!colnames(survivaldata) %in% "PATNUM"]

            survivaldata <- survivaldata[!survivaldata[,3] %in% c(NA, "") &
                                             !survivaldata[,2] %in% c(NA, "NotAvailable") & !survivaldata[,1] %in% c(NA, "NotAvailable"),]
            survivaldata[,1] <- as.numeric(survivaldata[,1])
            survivaldata[,2] <- as.numeric(survivaldata[,2])

            ## Ok, for the competing events - we have to clean those out cause it screws with the analysis - just remove I guess?
            survivaldata <- survivaldata[!survivaldata[,grepl("^C_", colnames(survivaldata))] %in% 2,]

            ## Need to properly factorize our cohort data
            if (COIlabel %in% c("IMGDEGIS")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("None", "Mild", "Moderate", "Severe"))
            } else if (COIlabel %in% c("IMGDEGIS_01_2_3")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Moderate", "Severe"))
            } else if (COIlabel %in% c("IMGDEGIS_01_3")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Severe"))
            } else if (COIlabel %in% c("CTNDV50")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("Non-evaluable", "1", "2", "3"))
            } else if (COIlabel %in% c("CTNDV50_123_clean")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("1", "2", "3"))
            } else if (COIlabel %in% c("CTNDV50_13_clean")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("1", "3"))
            } else if (COIlabel %in% c("CKD")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("No", "Yes"))
            } else if (COIlabel %in% c("highlowage")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("AGE_RAND_low", "AGE_RAND_highQ3"))
                
                
            } else if (COIlabel %in% c("rna_4cluster_w3AB")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3A", "nmf_cluster_3B", "nmf_cluster_4"))
                clusterplotCOI <- c(clusterplotCOI, "rna_4cluster_w3AB_nmf_cluster_1", "rna_4cluster_w3AB_nmf_cluster_2",
                                    "rna_4cluster_w3AB_nmf_cluster_3A", "rna_4cluster_w3AB_nmf_cluster_3B", "rna_4cluster_w3AB_nmf_cluster_4")
            } else if (COIlabel %in% c("rna_4cluster")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4"))
                clusterplotCOI <- c(clusterplotCOI, "rna_4cluster_nmf_cluster_1", "rna_4cluster_nmf_cluster_2",
                                    "rna_4cluster_nmf_cluster_3", "rna_4cluster_nmf_cluster_4")
            } else if (COIlabel %in% c("meth_4cluster")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("meth_4cluster_1", "meth_4cluster_2", "meth_4cluster_3", "meth_4cluster_4"))
                clusterplotCOI <- c(clusterplotCOI, "meth_4cluster_meth_4cluster_1", "meth_4cluster_meth_4cluster_2",
                                    "meth_4cluster_meth_4cluster_3", "meth_4cluster_meth_4cluster_4")
            } else if (COIlabel %in% c("meth_3cluster")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("meth_4cluster_1", "meth_4cluster_2", "meth_4cluster_3"))
                clusterplotCOI <- c(clusterplotCOI, "meth_4cluster_meth_4cluster_1", "meth_4cluster_meth_4cluster_2",
                                    "meth_4cluster_meth_4cluster_3")
            } else {
                cohortlabels <- unique(survivaldata[,3])
                survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]),
                                                                        as.character(cohortlabels[!grepl("not_", cohortlabels)])))
            }

            # At this point - lets also output an N table
            cohortNtablelist[[eventnum]] <- cbind(data.frame(table(survivaldata[,3])), COIlabel)
            names(cohortNtablelist)[eventnum] <- eventsel

            ## ADDING IN HERE A FAILSAFE - 2023-04-06
            if (length(table(survivaldata[,3])) < 2 | length(table(survivaldata[,3])) < 2) {
                outsurvtable <- data.frame(NA)
                outsurvplot <- NULL
                outsurvpvalue <- data.frame(variable = "Cohort", pval = 1, method = "Log-rank", pval.txt = 1)
            } else {
                ## Run the analysis
                survivalanalysisout <- suppressWarnings(create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL, ylimitparam = c(0.25,1)))
                outsurvtable <- survivalanalysisout$outsurvtable
                outsurvplot <- survivalanalysisout$outsurvplot
                outsurvpvalue <- survivalanalysisout$outsurvpvalue
            }
            ## ADDING IN HERE A FAILSAGE - 2023-04-06

            ## Save out the pvalue of the log-rank test
            survpvallist[[eventnum]] <- outsurvpvalue[,"pval"]
            names(survpvallist)[eventnum] <- eventsel

            # Ok, I think we need to return the HR as well for coxph to get some kind of effect size for this analysis
            outcoxphobject <- survivalanalysisout$outcoxphobject
            ## Put a hack in here to get the orientation correct. If theres only one comp, then coerce to a 1x3 dataframe:
            coxtempout <- summary(outcoxphobject)[["conf.int"]][,c("exp(coef)", "lower .95", "upper .95")]
            if(is.null(dim(coxtempout))){
                coxtempout <- t(data.frame(coxtempout))
                rownames(coxtempout) <- levels(survivaldata[,3])[2] ## I can do this because I know the first level is always the ref and I set this earlier above ^^
            }
            hrvalue <- cbind(coxtempout, Event = eventsel, Cohort = COIlabel)

            ## Save out the coxph HR of the coxph test
            survHRlist[[eventnum]] <- hrvalue
            names(survHRlist)[eventnum] <- eventsel

            ## Write out the plot and table
            outsubdir <- paste0(survivalplot_outfilepath, COIlabel, "/", eventsel, "/")
            dir.create(outsubdir, showWarnings = FALSE, recursive = TRUE)

            write.table(outsurvtable, paste0(outsubdir, COIlabel, "_", eventsel, "_survival_stat_table.csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            pdf(paste0(outsubdir, COIlabel, "_", eventsel, "_survival_plot.pdf"), width = 8, height = 5, onefile=FALSE)
            suppressWarnings(print(outsurvplot)) # warn suppress for when we cut off the stray point below 0.25
            junk <- dev.off()

        }
        survpvaltab <- do.call(rbind, survpvallist)
        colnames(survpvaltab) <- paste0(COIlabel, "_pval")
        survhazardtab <- do.call(rbind, survHRlist)
        colnames(survhazardtab) <- paste0(COIlabel, "_", colnames(survhazardtab))
        cohortNtable <- do.call(rbind, cohortNtablelist)
        colnames(cohortNtable) <- c("Cohort", "Number", "COI")


        fullpvaloutlist[[COInum]] <- survpvaltab
        names(fullpvaloutlist[COInum]) <- paste0(COIlabel)
        fullcumhazoutlist[[COInum]] <- survhazardtab
        names(fullcumhazoutlist)[COInum] <- paste0(COIlabel)
        fullcohortNtablelist[[COInum]] <- cohortNtable
        names(fullcohortNtablelist)[COInum] <- paste0(COIlabel)


    }
    fullpvalouttable <- do.call(cbind, fullpvaloutlist)
    write.table(fullpvalouttable, paste0(survivalplot_outfilepath, "full_survival_pval_table.csv"), sep = ",", col.names = NA, row.names = TRUE)

    fullcumhazouttable <- data.frame(do.call(rbind, fullcumhazoutlist))
    colnames(fullcumhazouttable) <- c("HR", "lower_CI", "upper_CI", "Event", "Cohort")
    fullcumhazouttable[,c("HR", "lower_CI", "upper_CI")] <- apply(fullcumhazouttable[,c("HR", "lower_CI", "upper_CI")], 2, function(x)
        as.numeric(as.character(x)))
    write.table(fullcumhazouttable, paste0(survivalplot_outfilepath, "full_survival_cumhaz_table.csv"), sep = ",", col.names = NA, row.names = TRUE)

    fullcohortNouttable <- do.call(rbind, fullcohortNtablelist)
    write.table(fullcohortNouttable, paste0(survivalplot_outfilepath, "full_survival_cohortN_table.csv"), sep = ",", col.names = NA, row.names = TRUE)

    ## Summary heatmap of the pvals
    fullstatouttable <- fullpvalouttable

    ## Turn the event outcome into a heatmap
    eventpvalcutoff <- eventpvalcutoff
    outeventpvaltab <- fullstatouttable[,paste0(plotCOI, "_pval")]

    ## Create the initial maptab
    maptab <- data.frame(outeventpvaltab)
    EOIvec <- colnames(outeventpvaltab)
    ## Add in geneordering
    geneorder <- rownames(outeventpvaltab[order(outeventpvaltab[,1], decreasing = FALSE),])

    ## Create HR table for directionality
    hazardmaptabtemp <- dcast(data.frame(fullcumhazouttable[fullcumhazouttable[,"Cohort"] %in% plotCOI, c("Event", "Cohort", "HR")]),
                              Event ~ Cohort, value.var = "HR")
    rownames(hazardmaptabtemp) <- hazardmaptabtemp[,1]
    hazardmaptabtemp <- hazardmaptabtemp[rownames(maptab), gsub("_pval", "", colnames(maptab))]
    hazardmaptab <- apply(hazardmaptabtemp, 2, as.numeric)
    rownames(hazardmaptab) <- hazardmaptabtemp[,1]
    hazardsigntab <- ifelse(hazardmaptab > 1, 1, -1)

    maptab <- maptab * hazardsigntab

    heatmapcolorparam = colorRamp2(c(-1, -eventpvalcutoff - 0.0001, -eventpvalcutoff, -0.0001, -1e-100, 0,
                                     1e-100, 0.0001, eventpvalcutoff, eventpvalcutoff + 0.0001, 1),
                                   c("grey", "grey", "#bcb2fd", "#11007d", "#11007d", "white",
                                     "#7d0000", "#7d0000", "#ffb2b2", "grey", "grey"))
    # "#bcb2fd"
    # old light red "#fd9999"
    # heatmapcolorLEGEND = colorRamp2(c(-eventpvalcutoff, 0, eventpvalcutoff), c("blue", "white", "red"))

    ht1 = Heatmap(as.matrix(maptab),
                  col = heatmapcolorparam,    ## Define the color scale for the heatmap
                  row_title = "Events", column_title = "Cohorts",
                  border = TRUE, na_col = "white", rect_gp = gpar(col = "black", lwd = 0.5),

                  cluster_columns = FALSE, cluster_rows = FALSE,clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                  #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"

                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize=6),
                  show_row_dend = TRUE, show_column_dend = TRUE,

                  heatmap_legend_param = list(
                      title = "KM pval", at = c(-0.05, 0, 0.05)),
                  height = unit(min((nrow(maptab)/2), 12),"cm"), width = unit(min(ncol(maptab), 18),"cm")

    )
    
    
    if (return_output_tables == TRUE) {
        return(list(plotout = ht1, tableout = maptab))
    } else  {
        return(ht1)
    }

}


# --------------------------------- cluster vs event - but nonbinary! --------------------------------- 

cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap <- function(selected_cluster_table, COIref_tablist,
                                                                  survivalplot_outfilepath, bioreptable_waddons, eventpvalcutoff, return_output_tables = FALSE,
                                                                  return_pairwise_coxph_values = FALSE, coxph_control_table = NULL, EOI = NULL, ...) {
    
    outfilepathsurvival <- survivalplot_outfilepath
    if(is.null(EOI)) {EOI <- gsub("^C_", "", colnames(bioreptable_waddons)[grepl("^C_", colnames(bioreptable_waddons))])}
    EOIcols <- apply(expand.grid(c("T_", "C_"), EOI, stringsAsFactors = FALSE), 1, paste, collapse = "")
    selected_cluster_table_list <- list(selected_cluster_table[,2,drop=FALSE])
    names(selected_cluster_table_list) <- colnames(selected_cluster_table[,2,drop=FALSE])
    
    ## Ok - so we want to iterate over every GOI, and every event, and see what comes out
    survival_fulltable <- bioreptable_waddons
    write.table(survival_fulltable, paste0(outfilepathsurvival, "full_survival_intable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    COItablist <- c(COIref_tablist, selected_cluster_table_list)
    fullpvaloutlist <- fullcumhazoutlist <- fullcohortNtablelist <- list()
    clusterplotCOI <- c()
    for (COInum in seq_len(length(COItablist))) {
        ## Select COI
        COIsel <- COItablist[[COInum]]
        COIlabel <- names(COItablist)[COInum]

        survival_intable <- merge(survival_fulltable[,c("PATNUM", EOIcols)], cbind(PATNUM = rownames(COIsel), COIsel), by = "PATNUM")
        colnames(survival_intable)[ncol(survival_intable)] <- COIlabel
        
        survpvallist <- survHRlist <- cohortNtablelist <- list()
        for (eventnum in seq_len(length(EOI))) {
            ## Select event
            eventsel <- EOI[eventnum]
            
            ## Create the subsurvivaltab
            # survivaldata <- data.frame(na.omit(survival_intable[,c(paste0(c("T_", "C_"), eventsel), COIlabel)]), row.names = survival_intable[,"PATNUM"])
            survivaldata <- na.omit(data.frame(survival_intable[,c(paste0(c("T_", "C_"), eventsel), COIlabel, "PATNUM")]))
            rownames(survivaldata) <- survivaldata[,"PATNUM"]
            survivaldata <- survivaldata[,!colnames(survivaldata) %in% "PATNUM"]
            
            survivaldata <- survivaldata[!survivaldata[,3] %in% c(NA, "") &
                                         !survivaldata[,2] %in% c(NA, "NotAvailable") & !survivaldata[,1] %in% c(NA, "NotAvailable"),]
            survivaldata[,1] <- as.numeric(survivaldata[,1])
            survivaldata[,2] <- as.numeric(survivaldata[,2])
            
            ## Ok, for the competing events - we have to clean those out cause it screws with the analysis - just remove I guess?
            survivaldata <- survivaldata[!survivaldata[,grepl("^C_", colnames(survivaldata))] %in% 2,]
            
            # ## Apply a max time of 5 years - it cleans up the data a little - NOPE, dont need it
            # maxtimelimit <- 5*365
            # survivaldata <- survivaldata[survivaldata[,grepl("^T_", colnames(survivaldata))] < maxtimelimit,]
            
            ## Need to properly factorize our cohort data
            if (COIlabel %in% c("IMGDEGIS")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("None", "Mild", "Moderate", "Severe"))
            } else if (COIlabel %in% c("IMGDEGIS_01_2_3")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Moderate", "Severe"))
            } else if (COIlabel %in% c("CUSTOM_IMGDEGIS_NONEMILD")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Moderate", "Severe"))
            } else if (COIlabel %in% c("IMGDEGIS_01_3")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Severe"))
            } else if (COIlabel %in% c("CTNDV50")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("Non-evaluable", "1", "2", "3"))
            } else if (COIlabel %in% c("CTNDV70")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("Non-evaluable", "0", "1", "2", "3"))
            } else if (COIlabel %in% c("CTNDV50_123_clean")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("1", "2", "3"))
            } else if (COIlabel %in% c("CTNDV50_13_clean")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("1", "3"))
            } else if (COIlabel %in% c("CKD")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("No", "Yes"))
            } else if (COIlabel %in% c("highlowage")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("AGE_RAND_low", "AGE_RAND_highQ3"))
            } else if (COIlabel %in% c("age_cat")) {
                survivaldata[,3] <- factor(survivaldata[,3], levels = c("<= 34", "35 - 44", "45 - 54", "55 - 64", "65 - 74", ">= 75"))
                
                
            } else if (COIlabel %in% c("rna_4cluster_w3AB")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3A", "nmf_cluster_3B", "nmf_cluster_4"))
                clusterplotCOI <- c(clusterplotCOI, "rna_4cluster_w3AB_nmf_cluster_1", "rna_4cluster_w3AB_nmf_cluster_2",
                                    "rna_4cluster_w3AB_nmf_cluster_3A", "rna_4cluster_w3AB_nmf_cluster_3B", "rna_4cluster_w3AB_nmf_cluster_4")
            } else if (COIlabel %in% c("rna_4cluster")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4"))
                clusterplotCOI <- c(clusterplotCOI, "rna_4cluster_nmf_cluster_1", "rna_4cluster_nmf_cluster_2",
                                    "rna_4cluster_nmf_cluster_3", "rna_4cluster_nmf_cluster_4")
            } else if (COIlabel %in% c("meth_4cluster")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("meth_4cluster_1", "meth_4cluster_2", "meth_4cluster_3", "meth_4cluster_4"))
                clusterplotCOI <- c(clusterplotCOI, "meth_4cluster_meth_4cluster_1", "meth_4cluster_meth_4cluster_2",
                                    "meth_4cluster_meth_4cluster_3", "meth_4cluster_meth_4cluster_4")
            } else if (COIlabel %in% c("meth_3cluster")) {
                survivaldata[,3] <- factor(survivaldata[,3], 
                                           levels = c("meth_3cluster_1", "meth_3cluster_2", "meth_3cluster_3"))
                clusterplotCOI <- c(clusterplotCOI, "meth_3cluster_meth_3cluster_1", "meth_3cluster_meth_3cluster_2",
                                    "meth_3cluster_meth_3cluster_3")
            } else {
                cohortlabels <- unique(survivaldata[,3])
                survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]),
                                                                        as.character(cohortlabels[!grepl("not_", cohortlabels)])))
            }
            
            # At this point - lets also output an N table
            cohortNtablelist[[eventnum]] <- cbind(data.frame(table(survivaldata[,3])), COIlabel)
            names(cohortNtablelist)[eventnum] <- eventsel
            
            # ## Run the analysis
            # survivalanalysisout <- suppressWarnings(create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL, ylimitparam = c(0.25,1)))
            # outsurvtable <- survivalanalysisout$outsurvtable
            # outsurvplot <- survivalanalysisout$outsurvplot
            # outsurvpvalue <- survivalanalysisout$outsurvpvalue
            
            ## ADDING IN HERE A FAILSAFE - 2023-04-06
            if (length(table(survivaldata[,3])) < 2 | length(table(survivaldata[,3])) < 2) {
                outsurvtable <- data.frame(NA)
                outsurvplot <- NULL
                outsurvpvalue <- data.frame(variable = "Cohort", pval = 1, method = "Log-rank", pval.txt = 1)
            } else {
                ## Run the analysis
                survivalanalysisout <- suppressWarnings(create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL, ylimitparam = c(0.25,1)))
                outsurvtable <- survivalanalysisout$outsurvtable
                outsurvplot <- survivalanalysisout$outsurvplot
                outsurvpvalue <- survivalanalysisout$outsurvpvalue
            }
            ## ADDING IN HERE A FAILSAGE - 2023-04-06
            
            ## Save out the pvalue of the log-rank test
            survpvallist[[eventnum]] <- outsurvpvalue[,"pval"]
            names(survpvallist)[eventnum] <- eventsel
            
            # So here - I can grab the pairwise Pvals as well if I want..., but is that the same as the log-rank of a pairwise? Maybe, maybe not
            ## Actually - this also only returns the first factor vs all of the non first factor, and not the OTHER ovals as well...
            if (return_pairwise_coxph_values) {
                outsubdir <- paste0(outfilepathsurvival, COIlabel, "/", eventsel, "/")
                dir.create(outsubdir, showWarnings = FALSE, recursive = TRUE)
                
                modeltable <- survivaldata
                colnames(modeltable) <- c("Time_to_event", "event_status", "Cohort")
                class_labels <- unique(modeltable[,3])
                pairwise_coxout_list <- list()
                for (class_label in class_labels) {
                    temp_modeltable <- modeltable
                    temp_modeltable[,3] <- factor(temp_modeltable[,3], levels = c(as.factor(class_label), 
                                                                                  sort(class_labels[!class_labels %in% class_label])))
                    
                    if (!is.null(coxph_control_table)) {
                        temp_modeltable <- cbind(temp_modeltable, coxph_control_table[rownames(temp_modeltable),])
                        temp_coxobject <- suppressWarnings(coxph(as.formula(paste0("Surv(Time_to_event, event_status) ~ ", 
                                                                  paste(colnames(temp_modeltable)[3:(ncol(temp_modeltable))], collapse = " + "))), 
                                                data = temp_modeltable))
                    } else {

                        temp_coxobject <- coxph(Surv(Time_to_event, event_status) ~ Cohort, data = temp_modeltable)
                    }
                    
                    out1 <- cbind.data.frame(summary(temp_coxobject)[["coefficients"]], summary(temp_coxobject)[["conf.int"]])
                    if ("xlevels" %in% names(temp_coxobject)) {
                        out1[,"control_class"] <- unlist(lapply(temp_coxobject[["xlevels"]], function(x) rep(x[1], (length(x)-1))))
                    } else {
                        out1[,"control_class"] <- class_label
                    }
                    
                    out1[,"test_class"] <- rownames(out1)
                    out1[,"iteration"] <- class_label
                    pairwise_coxout_list[[class_label]] <- out1
                }
                pairwise_coxout_table <- do.call(rbind, pairwise_coxout_list)
                
                # Now just grab the values that we actually care about
                coxph_pair_pval_table <- reshape2::dcast(pairwise_coxout_table[pairwise_coxout_table[,"control_class"] %in% class_labels,], 
                                                         formula = control_class ~ test_class, value.var = "Pr(>|z|)")
                plotout_table <- data.frame(coxph_pair_pval_table[,2:ncol(coxph_pair_pval_table)], row.names = coxph_pair_pval_table[,1])
                plotout_table <- round(plotout_table, 3)
                
                write.table(pairwise_coxout_table, paste0(outsubdir, COIlabel, "_", eventsel, "_pairwise_coxph_table.csv"),
                            sep = ",", col.names = TRUE, row.names = FALSE)
                pdf(paste0(outsubdir, COIlabel, "_", eventsel, "_pairwise_coxph_tableplot.pdf"), width = 10, height = 5, onefile=FALSE)
                plot(tableGrob(plotout_table,
                               theme = ttheme(base_style = "classic", 
                                              rownames.style = rownames_style(fontface = "bold"))
                               )
                     )
                junk <- dev.off()
            }

            
            # I think there is a way to do this otherwise, but for now - lets just run coxph with each factor as baseline, then concat and clean
            
            
            # Ok, I think we need to return the HR as well for coxph to get some kind of effect size for this analysis
            outcoxphobject <- survivalanalysisout$outcoxphobject
            ## Put a hack in here to get the orientation correct. If theres only one comp, then coerce to a 1x3 dataframe:
            coxtempout <- summary(outcoxphobject)[["conf.int"]][,c("exp(coef)", "lower .95", "upper .95")]
            if(is.null(dim(coxtempout))){
                coxtempout <- t(data.frame(coxtempout))
                rownames(coxtempout) <- levels(survivaldata[,3])[2] ## I can do this because I know the first level is always the ref and I set this earlier above ^^
            }
            hrvalue <- cbind(coxtempout, Event = eventsel, Cohort = COIlabel)
            
            ## Save out the coxph HR of the coxph test
            survHRlist[[eventnum]] <- hrvalue
            names(survHRlist)[eventnum] <- eventsel
            
            ## Write out the plot and table
            outsubdir <- paste0(outfilepathsurvival, COIlabel, "/", eventsel, "/")
            dir.create(outsubdir, showWarnings = FALSE, recursive = TRUE)
            
            write.table(outsurvtable, paste0(outsubdir, COIlabel, "_", eventsel, "_survival_stat_table.csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            pdf(paste0(outsubdir, COIlabel, "_", eventsel, "_survival_plot.pdf"), width = 8, height = 5, onefile=FALSE)
            suppressWarnings(print(outsurvplot))
            junk <- dev.off()
            
        }
        survpvaltab <- do.call(rbind, survpvallist)
        colnames(survpvaltab) <- paste0(COIlabel, "_pval")
        survhazardtab <- do.call(rbind, survHRlist)
        colnames(survhazardtab) <- paste0(COIlabel, "_", colnames(survhazardtab))
        cohortNtable <- do.call(rbind, cohortNtablelist)
        colnames(cohortNtable) <- c("Cohort", "Number", "COI")
        
        
        fullpvaloutlist[[COInum]] <- survpvaltab
        names(fullpvaloutlist[COInum]) <- paste0(COIlabel)
        fullcumhazoutlist[[COInum]] <- survhazardtab
        names(fullcumhazoutlist)[COInum] <- paste0(COIlabel)
        fullcohortNtablelist[[COInum]] <- cohortNtable
        names(fullcohortNtablelist)[COInum] <- paste0(COIlabel)
        
        
    }
    fullpvalouttable <- do.call(cbind, fullpvaloutlist)
    write.table(fullpvalouttable, paste0(outfilepathsurvival, "full_survival_pval_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    fullcumhazouttable <- data.frame(do.call(rbind, fullcumhazoutlist))
    colnames(fullcumhazouttable) <- c("HR", "lower_CI", "upper_CI", "Event", "Cohort")
    fullcumhazouttable[,c("HR", "lower_CI", "upper_CI")] <- apply(fullcumhazouttable[,c("HR", "lower_CI", "upper_CI")], 2, function(x)
        as.numeric(as.character(x)))
    write.table(fullcumhazouttable, paste0(outfilepathsurvival, "full_survival_cumhaz_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    fullcohortNouttable <- do.call(rbind, fullcohortNtablelist)
    write.table(fullcohortNouttable, paste0(outfilepathsurvival, "full_survival_cohortN_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    ## Summary heatmap of the pvals
    # fullstatouttablefile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run1c_rmoutliers2/integrative_analyses/nmf_survival_analysis/full_survival_pval_table.csv"
    # fullstatouttable <- read.table(fullstatouttablefile, sep = ",", header = TRUE, row.names = 1)
    fullstatouttable <- fullpvalouttable
    
    ## Turn the event outcome into a heatmap
    outeventpvaltab <- fullstatouttable[,paste0(plotCOI, "_pval")]
    maptab <- data.frame(outeventpvaltab)
    
    ## MULTIPLE HYPOTHESIS CORRECTION - OMITTING FOR NOW BECAUSE WE HAVE WAY TOO MANY BAD HYPOTHESES AND ITS CLOGGING IT UP
    # p.adjust(unlist(maptab), method = "fdr") ## Sanity check - it works
    # “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
    
    # correctedmat <- matrix(p.adjust(as.vector(as.matrix(maptab)), method='fdr'),ncol=ncol(maptab))
    # rownames(correctedmat) <- rownames(maptab)
    # colnames(correctedmat) <- colnames(maptab)
    # maptab <- correctedmat
    
    # Create a heatmap of just the event pvalues
    EOIvec <- colnames(outeventpvaltab)
    # eventorder <- rownames(outeventpvaltab[order(outeventpvaltab[,1], decreasing = FALSE),])
    heatmapcolorparam = colorRamp2(c(-1, -eventpvalcutoff - 0.0001, -eventpvalcutoff, -0.0001, -1e-100, 0,
                                     1e-100, 0.0001, eventpvalcutoff, eventpvalcutoff + 0.0001, 1),
                                   c("grey", "grey", "#bcb2fd", "#11007d", "#11007d", "white",
                                     "#7d0000", "#7d0000", "#ffb2b2", "grey", "grey"))
    
    ht1 = Heatmap(as.matrix(maptab),
                  col = heatmapcolorparam,    ## Define the color scale for the heatmap
                  row_title = "Events", column_title = "Cohorts",
                  border = TRUE, na_col = "white", rect_gp = gpar(col = "black", lwd = 0.5),
                  
                  cluster_columns = FALSE, cluster_rows = FALSE,clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                  #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
                  
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize=6),
                  show_row_dend = TRUE, show_column_dend = TRUE,
                  
                  # heatmap_legend_param = list(
                  #     col_fun = heatmapcolorLEGEND,
                  #     title = "pvalue", legend_height = unit(2.5, "cm"), title_gp = gpar(fontsize = 8, fontface = "bold")),
                  heatmap_legend_param = list(
                      title = "KM pval", at = c(-0.05, 0, 0.05)),
                  height = unit(min((nrow(maptab)/2), 12),"cm"), width = unit(min(ncol(maptab), 18),"cm")
                  
    )
    
    
    if (return_output_tables == TRUE) {
        return(list(plotout = ht1, tableout = maptab))
    } else  {
        return(ht1)
    }
    
}




# --------------------------------- subtype deseq summary heatmap --------------------------------- 

# parameter added for which deseq genes: all, positive, or negative
subtype_deseq_summary_heatmap <- function(selected_deseq_table_list, whichgenes = "all",
                                          normcounttable, selected_rna_clustermembership_table, 
                                          pval_test = "padj", pval_cutoff = 0.05, log2fc_cutoff = 0.5,
                                          rowclusterparam = TRUE, ...) {
    gene_annotation_table_list <- list()
    for (subtype in names(selected_deseq_table_list)) {
        table_select <- selected_deseq_table_list[[subtype]]
        name_select <- strsplit(subtype, split = "_")[[1]][1]
        outtable <- rbind(cbind.data.frame(genes = rownames(table_select[table_select[, pval_test] < pval_cutoff & 
                                                                         !is.na(table_select[, pval_test]) &
                                                                         table_select[,"log2FoldChange"] > log2fc_cutoff,]),
                                           subtype = 1),
                          cbind.data.frame(genes = rownames(table_select[table_select[, pval_test] < pval_cutoff &
                                                                         !is.na(table_select[, pval_test]) &
                                                                         table_select[,"log2FoldChange"] < -log2fc_cutoff,]),
                                           subtype = -1))
        colnames(outtable)[2] <- name_select
        gene_annotation_table_list[[name_select]] <- outtable
    }
    # I think the issue is the 0s and 1s and -1s, and we need JUST 1s and 0s
    # apply(gene_annotation_table[,2:ncol(gene_annotation_table)], 2, function(x) table(x, useNA = "always"))
    
    gene_annotation_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "genes", all = TRUE, sort = FALSE), gene_annotation_table_list)
    rownames(gene_annotation_table) <- gene_annotation_table[,"genes"]
    gene_annotation_table <- gene_annotation_table[,!grepl("genes", colnames(gene_annotation_table))]
    
    # Really I just need all of this conversion nonsense for the geneorder and thats it
    # Add in 0s, turn all we dont care about to 0, sort, and then keep that gene order
    gene_order_table <- gene_annotation_table
    gene_order_table[is.na(gene_order_table)] <- 0
    # apply(gene_order_table[,2:ncol(gene_order_table)], 2, function(x) table(x, useNA = "always"))
    if (whichgenes == "positive") {
        gene_order_table[gene_order_table != 1] <- 0
    } else if (whichgenes == "negative") {
        gene_order_table[gene_order_table != -1] <- 0
        gene_order_table[gene_order_table == -1] <- 1
    }
    gene_order_table <- gene_order_table[order(-gene_order_table[,1], -gene_order_table[,2],
                                               -gene_order_table[,3], -gene_order_table[,4]), ]
    GOI <- rownames(gene_order_table)

    
    # Turn back into text
    gene_annotation_table[gene_annotation_table == 1] <- "positive"
    gene_annotation_table[gene_annotation_table == -1] <- "negative"
    gene_annotation_table[gene_annotation_table == 0] <- NA
    
    deseq_summary_hm_plottab <- normcounttable[GOI, colnames(normcounttable)]
    colmetatable <- selected_rna_clustermembership_table[colnames(normcounttable),,drop=FALSE]
    colannotationlist <- annotationlist_builder(colmetatable)
    rowmetatable <- gene_annotation_table[GOI,]
    
    customcolorlist <- rep(list(c(positive = "#E31A1C", negative = "dodgerblue2")), ncol(rowmetatable))
    names(customcolorlist) <- colnames(rowmetatable)
    rowannotationlist <- annotationlist_builder(rowmetatable, customcolorlist = customcolorlist)
    
    # Which genes?
    if (whichgenes == "positive") {
        selected_genes <- rownames(rowmetatable[rowSums(rowmetatable == "positive", na.rm = TRUE) > 0,])
        rowmetatable <- rowmetatable[GOI[GOI %in% selected_genes],]
        rowmetatable[rowmetatable == "negative"] <- NA
    } else if (whichgenes == "negative") {
        selected_genes <- rownames(rowmetatable[rowSums(rowmetatable == "negative", na.rm = TRUE) > 0,])
        rowmetatable <- rowmetatable[GOI[GOI %in% selected_genes],]
        rowmetatable[rowmetatable == "positive"] <- NA
    } else {
        selected_genes <- rownames(rowmetatable)
    }
    
    
    # Manually scaling to shrink the range
    heatmapcolorparam = colorRamp2(breaks = c(-5, 0, 5), colors = c("blue", "white", "red"))
    
    outhm1 <- create_heatmap(counttab = deseq_summary_hm_plottab[GOI[GOI %in% selected_genes],], scale_data = TRUE, separate_legend = FALSE,
                             colmetatable = colmetatable, colannotationlist = colannotationlist,
                             rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                             colclusterparam = TRUE, rowclusterparam = rowclusterparam,
                             columnsplitparam = colmetatable,
                             heatmapcolorparam = heatmapcolorparam)
    return(outhm1$heatmap)
    # draw(outhm1$heatmap)
    # outheatmapfile1 <- paste0(outfilepathmaster, "nmfcluster_alldiffexpgenes_heatmap.pdf")
    # pdf(outheatmapfile1, useDingbats = FALSE, height = 7, width = 13)
    # draw(outhm1$heatmap)
    # junk <- dev.off()
}

# --------------------------------- clean volcano plot function --------------------------------- 

clean_volcano_plot_function <- function(deseqtable, nameparam="volcano", labeledgenes = NULL, 
                                        pval_test = "padj", pval_cutoff = 0.05, log2fc_cutoff = 0.5, returnplottab = FALSE) {
    stattype = pval_test
    pvalcutoff = pval_cutoff
    log2fccutoff = log2fc_cutoff
    
    ## PLOT CODE
    plottab <- deseqtable[,c("log2FoldChange", stattype)]
    plottab[is.na(plottab[,2]),2] <- 1 ## Failsafe to have any pvalues of NA be set as 1
    plottab[is.na(plottab[,1]),1] <- 0
    
    plottab$color = "notsig"
    plottab[plottab[,2] < pvalcutoff & plottab[,1,drop=FALSE] > 0, "color"] <- "sigpvalup" # pink
    plottab[plottab[,1,drop=FALSE] > log2fccutoff, "color"] <- "siglog2fcup" # red
    plottab[plottab[,1,drop=FALSE] > log2fccutoff & plottab[,2] < pvalcutoff, "color"] <- "sigbothup" # darkred
    
    plottab[plottab[,2] < pvalcutoff & plottab[,1,drop=FALSE] < 0, "color"] <- "sigpvaldown" # light blue
    plottab[plottab[,1,drop=FALSE] < -log2fccutoff, "color"] <- "siglog2fcdown" # blue
    plottab[plottab[,1,drop=FALSE] < -log2fccutoff & plottab[,2] < pvalcutoff, "color"] <- "sigbothdown" # dark blue
    
    plottab[,2] = -log10(plottab[,2])
    plottabin <- as.data.frame(plottab)
    # great, pink, red, darkred, light blue, blue, dark blue
    # colorparam = c("#e5e5e5", "#ff9999", "#ff0000", "#990000", "#b2b2ff", "#1919ff", "#0000b2") 
    colorparam = c("#e5e5e5", "#ff9999", "#ff9999", "#990000", "#b2b2ff", "#b2b2ff", "#0000b2") 
    names(colorparam) = c("notsig","sigpvalup", "siglog2fcup", "sigbothup", "sigpvaldown", "siglog2fcdown", "sigbothdown" )
    
    # # Original version - full pdf
    # pout <- ggplot(plottabin, aes(x = plottabin[,1], y = plottabin[,2], color = plottabin[,3]))
    # pout <- pout + geom_point(size = 1) + scale_color_manual(values = colorparam)
    # pout <- pout + labs(x = "Log2FoldChange", y = "-Log10(P. Value)", title = nameparam, color = "significance") + theme_pubr(base_size = 20)
    # pout <- pout + theme(legend.position = "none")
    # 
    # 
    # labeltab = plottabin[rownames(plottabin) %in% labeledgenes,]
    # pout <- pout + geom_label_repel(max.overlaps = 100,
    #                                 data = labeltab,
    #                                 aes(x = labeltab[,1], y=labeltab[,2], label = rownames(labeltab), color=labeltab[,3]),
    #                                 size=8, segment.size = 0.2, fontface="bold", segment.color = "black",
    #                                 box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"))
    
                                    
    # Rasterized version for BIG volcano plots
    # Set colors first before running:
    pout <- ggplot(plottabin, aes(x = plottabin[,1], y = plottabin[,2]))
    temp_color_table <- plottabin[,3,drop=FALSE]
    colnames(temp_color_table) <- "colorvar"
    temp_color_table[,"keeporder"] <- seq(1, nrow(temp_color_table))
    tt1 <- merge(temp_color_table, data.frame(colorparam),
                 by.x = "colorvar", by.y = "row.names", sort = FALSE)
    tt2 <- tt1[order(tt1[,"keeporder"]),]
    color_for_plotting <- tt2[,3,drop = FALSE]
    colnames(color_for_plotting) <- "colorvar"

    pout <- pout + geom_scattermore(color = color_for_plotting[,1], fill = color_for_plotting[,1], pointsize = 3, pixels = c(700,700), shape = 19)

    # labeltab <- merge(cbind(feature_label = rownames(plottabin[rownames(plottabin) %in% labeledgenes,]),
    #                         plottabin[rownames(plottabin) %in% labeledgenes,]),
    #                   data.frame(colorparam), by.x = "color", by.y = "row.names", sort = FALSE)
    labeltab = plottabin[rownames(plottabin) %in% labeledgenes,]
    pout <- pout + geom_label_repel(max.overlaps = 100,
                                    data = labeltab,
                                    aes(x = labeltab[,"log2FoldChange"], y=labeltab[,stattype], label = rownames(labeltab), color=labeltab[,"color"]),
                                    size = 12, segment.size = 0.2, fontface="bold", segment.color = "black",
                                    box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"))
    pout <- pout + scale_colour_manual(values = colorparam)
    
    pout <- pout + labs(x = "Log2FoldChange", y = "-Log10(P. Value)", title = nameparam, color = "significance") + theme_pubr(base_size = 20)
    pout <- pout + theme(legend.position = "none")
    pout
    
    if (returnplottab) {
        return(list(pout = pout, plottabin = plottabin))
    } else {
        return(pout)
    }
    
}


# --------------------------------- clean gsea barplot function --------------------------------- 

clean_gsea_plot_function <- function(gseatable, nameparam, pathwayselect,
                                     pval_test = "p.adjust", pval_cutoff = 0.05) {
    stattype = pval_test
    pvalcutoff = pval_cutoff

    POI <- pathwayselect
    
    titleparam <- strsplit(nameparam, split = "_")[[1]][1]
    if (nameparam == "RNAtype3B_vs_RNAtype3A") {titleparam <- "RNAtype3Bv3A"}
    
    plottab <- gseatable[POI, c("ID", "NES", stattype)]
    plottab <- plottab[order(plottab[,"NES"], decreasing = FALSE),]
    plottab[,1] <- gsub("GOMF_|GOBP_|GOCC_", "", plottab[,1])
    plottab[,1] <- sapply(tolower(gsub("_", " ", plottab[,1])), simpleCap)
    plottab[,"colorvar"] <- ifelse(plottab[,"NES"] > 0, "darkred", "darkblue")
    
    GOtermpout <- plot_barchart(bartab = plottab[,c(1,2)], colorvar = plottab[,"colorvar",drop=FALSE], 
                                labsparam = list(title = titleparam, x = "GO Term", y = "NES"),
                                strwrap_characternum = 1000)
    GOtermpout <- GOtermpout + coord_flip()
    GOtermpout <- GOtermpout + theme_pubr(base_size = 20)
    return(GOtermpout)
}




create_color_reference <- function(COI, colorguide) {
    COIlist <- list()
    for (COInum in seq_len(length(COI))) {
        COIsel <- COI[COInum]
        COIlist[[COInum]] <- list(ORDER = colorguide[colorguide[,"category"] %in% COIsel,"feature"],
                                  COLOR = colorguide[colorguide[,"category"] %in% COIsel,"color"])
        names(COIlist)[COInum] <- COIsel
    }
    return(COIlist)
}

# --------------------------------- rectangel plot for metadata features --------------------------------- 

create_rect_plot <- function(bioreptable_sel, COI_colref_list=NULL, plotaspercents = FALSE, addlabels = FALSE, addnumbers = FALSE) {
    subtable <- table(bioreptable_sel, useNA = "ifany")
    
    rectangleplot_table <- data.frame(subtable)
    rectangleplot_table[,"start"] <- c(0,cumsum(rectangleplot_table[1:(nrow(rectangleplot_table)-1),"Freq"]))
    rectangleplot_table[,"end"] <- cumsum(rectangleplot_table[1:(nrow(rectangleplot_table)),"Freq"])
    
    # Add more info so we can add text
    rectangleplot_table[,"mid"] <- (rectangleplot_table[,"end"] - rectangleplot_table[,"start"])/2 + rectangleplot_table[,"start"]
    
    # Make the rectangle plot out of 100% instead?
    if (plotaspercents == TRUE) {
        rectangleplot_table[,2:ncol(rectangleplot_table)] <- round(rectangleplot_table[,2:ncol(rectangleplot_table)] / sum(subtable), 3) * 100
    }
    
    rectplot_out <- ggplot(rectangleplot_table, aes(xmin = start, xmax = end, ymin = 0, ymax = 1))
    rectplot_out <- rectplot_out + geom_rect(aes(fill = bioreptable_sel), colour = "grey50")
    if (!is.null(COI_colref_list)) {
        rectplot_out <- rectplot_out + scale_fill_manual(values = COI_colref_list[[COI_sel]][["COLOR"]])
    }
    rectplot_out <- rectplot_out + theme_classic()
    
    # Add labels and numbers
    if (addlabels) {
        rectplot_out <- rectplot_out + geom_text(aes(x = mid, y = 0.25, label = bioreptable_sel))
    }
    if (addnumbers) {
        rectplot_out <- rectplot_out + geom_text(aes(x = mid, y = 0.75, label = Freq))
    }
    return(rectplot_out)
    
}

# --------------------------------- object to help clean up the DMP outputs (for the subtype outputs) --------------------------------- 
# TOSH YOU WANT DELTA BETA - notlof2fc it really doest make sense

clean_meth_subtype_DMP_object_function <- function(meth_subtype_DMP_object) {
    # DMP_labels <- names(meth_subtype_DMP_object)
    outlist <- list()
    for (DMP_table_num in seq_len(length(meth_subtype_DMP_object))) {
        grab_DMP_table <- meth_subtype_DMP_object[[DMP_table_num]][[1]]
        grab_label <- names(meth_subtype_DMP_object[[DMP_table_num]])
        # The logFC is a mess - so im just gonna totally redo and make it log2fc while im on it
        # theres the type_AVG and the not_AVG for all of these, and i want the type / not ratio
        notlabel <- "not_AVG"
        subtypelabel <- colnames(grab_DMP_table)[grepl("_AVG", colnames(grab_DMP_table)) & !grepl("not_", colnames(grab_DMP_table))]
        grab_DMP_table[,"signed_deltaBeta"] <- grab_DMP_table[,subtypelabel] - grab_DMP_table[,notlabel]
        # grab_DMP_table <- rename(grab_DMP_table, log2FoldChange = logFC)
        outlist[[DMP_table_num]] <- grab_DMP_table
        names(outlist)[DMP_table_num] <- grab_label
    }
    return(outlist)
}


# --------------------------------- write out a confusion matrix - helper function --------------------------------- 
## REPLACING OUT TEST AND TRAIN WITH FULL? TO GET BETTER CURVES
## OR MAYBE JUST WITH OUR TRIANING?
writeout_confmat <- function(confmat) {
    val1 <- data.frame(val = confmat$overall)
    val2 <- data.frame(val = confmat$byClass)
    val3 <- data.frame(val = confmat$table)
    rownames(val3) <- paste0("pred_", val3[,1], "_obs_", val3[,2])
    val3b <- val3[,3,drop=FALSE]
    colnames(val3b) <- "val"
    outtab <- do.call(rbind, list(val1, val2, val3b))
    return(outtab)
}

# --------------------------------- cohort classification ML pipeline --------------------------------- 
cohort_classification_ml_analysis <- function(featuretable, outcometable, outcomelevels, seedparam=11111,
                                              presplitdata = NULL,
                                              subsampleparam = NULL, models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")) {
    ## ML for SLE using Genes
    # packagelist = c("mlbench", "caret", "pROC")
    # junk <- lapply(packagelist, function(xxx) suppressMessages(
    #   require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
    # dir.create(paste0(modeling_outfilepath, "modeling/"), showWarnings = FALSE, recursive = TRUE)
    
    # Create modeling input data
    if (is.null(presplitdata)) {
        outcomelabel <- colnames(outcometable)
        featurelabels <- colnames(featuretable)
        modeling_intable <- na.omit(merge(outcometable, featuretable, by = "row.names"))
        rownames(modeling_intable) <- modeling_intable[,"Row.names"]
        modeling_intable <- modeling_intable[,!grepl("Row.names", colnames(modeling_intable))]
        modeling_intable[,outcomelabel] <- factor(modeling_intable[,outcomelabel], levels = outcomelevels)
        
        # Ok - we need to adjust variable names where applicable (like with 3v and 1v - doesnt like numbers first)
        level_key <- setNames(make.names(make.names(outcomelevels)), outcomelevels)
        modeling_intable[,outcomelabel] <- recode(modeling_intable[,outcomelabel], !!!level_key)
        outcomelevels <- make.names(outcomelevels)
        
        
        # create data partition and prepare training scheme
        #1415 0.92
        set.seed(seedparam)
        trainindex <- createDataPartition(modeling_intable[,outcomelabel], p = 0.6, list = FALSE)
        training <- modeling_intable[trainindex,]
        testing <- modeling_intable[-trainindex,]
    } else {
        training <- presplitdata[["training"]]
        testing <- presplitdata[["testing"]]
        outcomelabel <- colnames(training)[1]
        featurelabels <- colnames(training)[2:ncol(training)]
        
        # Still have to do the name catch/change - for training and testing separately
        level_key <- setNames(make.names(make.names(outcomelevels)), outcomelevels)
        training[,outcomelabel] <- recode(training[,outcomelabel], !!!level_key)
        testing[,outcomelabel] <- recode(testing[,outcomelabel], !!!level_key)
        outcomelevels <- make.names(outcomelevels)
        training[,outcomelabel] <- factor(training[,outcomelabel], levels = outcomelevels)
        testing[,outcomelabel] <- factor(testing[,outcomelabel], levels = outcomelevels)
    }
    # trainindex <- createDataPartition(modeling_intable[,outcomelabel], p = 1, list = FALSE)
    # training <- modeling_intable[trainindex,]
    # testing <- modeling_intable[trainindex,]
    if (is.null(subsampleparam)) {
        trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, 
                                  classProbs = TRUE, savePredictions = TRUE)
    } else {
        trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, 
                              classProbs = TRUE, savePredictions = TRUE,
                              sampling = subsampleparam  ## adding downsampling here
                              )
    }

    
    cohort_class_modeltraining <- function(training, testing, outcomelabel, featurelabels, trcontrol, modellabel, seedparam) {
        set.seed(seedparam)
        print(paste0("Training ", modellabel, " model"))

        if (modellabel == "xgbTree") {  ## Special params for xgboost
            modfit <- train(x=training[,featurelabels, drop=FALSE], y = training[,outcomelabel], method = modellabel, 
                            # trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC", verbose = FALSE, verbosity = 0)
                            trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC", verbose = FALSE, verbosity = 0)
        } else if (modellabel == "glmnet") { ## Special trick for lasso with a single variable
            training[,"ones"] <- rep(1, nrow(training))
            testing[,"ones"] <- rep(1, nrow(testing))
            featurelabels <- c(featurelabels, "ones")
            
            # A bunch of custom lambda coding:
            ## https://stackoverflow.com/questions/48280074/r-how-to-let-glmnet-select-lambda-while-providing-an-alpha-range-in-caret
            lambda <- unique(unlist(lapply(seq(0, 1, 0.05), function(x){
                init <- glmnet::glmnet(Matrix::as.matrix(training[,c(featurelabels)]), training[,outcomelabel],
                                       family = "binomial",
                                       nlambda = 100,
                                       alpha = x)
                lambda <- c(min(init$lambda), max(init$lambda))
            }
            )))
            tuneGridb <- expand.grid(.alpha = seq(0, 1, 0.05),
                                     .lambda = seq(min(lambda), max(lambda), length.out = 5))
            
            
            modfit <- suppressWarnings(train(x=training[,c(featurelabels), drop=FALSE], y = training[,outcomelabel], method = modellabel, 
                            trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC", verbose = FALSE, verbosity = 0, tuneGrid = tuneGridb))
        } else {
            modfit <- train(x=training[,featurelabels, drop=FALSE], y = training[,outcomelabel], method = modellabel, 
                            trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC")
        }
        predict_labels <- predict(modfit, testing[,featurelabels, drop=FALSE])
        predict_values <- predict(modfit, testing[,featurelabels, drop=FALSE], type = "prob")[,1,drop=FALSE]
        confstatsout <- writeout_confmat(confusionMatrix(predict_labels, testing[,outcomelabel]))
        predicted_values_and_labels <- cbind(predict_values, predict_labels)
        colnames(predicted_values_and_labels) <- c(paste0("predict_", modellabel, "_values"), paste0("predict_", modellabel, "_labels"))
        return(list(confstatsout = confstatsout, predicted_values_and_labels = predicted_values_and_labels))
    }
    # train some models
    
    # Im gonna do them separately instead of a for loop for now - probably not as neat, but easier to track for now
    models_ran <- c()
    if ("glm" %in% models_to_run) {
        glm_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                              modellabel = "glm", seedparam = seedparam)
        models_ran <- c(models_ran, "glm")
    } else { glm_out <- NULL }
    if ("xgbTree" %in% models_to_run) {
        xgbTree_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                  modellabel = "xgbTree", seedparam = seedparam) ## xgboost
        models_ran <- c(models_ran, "xgbTree")
    } else { xgbTree_out <- NULL }
    if ("svmRadial" %in% models_to_run) {
        svmRadial_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                    modellabel = "svmRadial", seedparam = seedparam)
        models_ran <- c(models_ran, "svmRadial")
    } else { svmRadial_out <- NULL }
    if ("rf" %in% models_to_run) {
        rf_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                             modellabel = "rf", seedparam = seedparam)
        models_ran <- c(models_ran, "rf")
    } else { rf_out <- NULL }
    if ("glmnet" %in% models_to_run) {
        glmnet_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                 modellabel = "glmnet", seedparam = seedparam) ## lasso
        models_ran <- c(models_ran, "glmnet")
    } else { glmnet_out <- NULL }
    if ("nnet" %in% models_to_run) {
        nnet_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                 modellabel = "nnet", seedparam = seedparam) ## neural net
        models_ran <- c(models_ran, "nnet")
    } else { nnet_out <- NULL }

    
    ## Attach to inmetatable by merge (for safety) then reassign and remove rownames
    modelstats_outtable = Reduce(function(dtf1, dtf2)
        merge(dtf1, dtf2, by = "metric", all.x = TRUE, sort = FALSE),
        Filter(Negate(is.null), list(
            cbind(metric = rownames(glm_out$confstatsout), glm_out$confstatsout),
            cbind(metric = rownames(xgbTree_out$confstatsout), xgbTree_out$confstatsout),
            cbind(metric = rownames(svmRadial_out$confstatsout), svmRadial_out$confstatsout), 
            cbind(metric = rownames(rf_out$confstatsout), rf_out$confstatsout), 
            cbind(metric = rownames(glmnet_out$confstatsout), glmnet_out$confstatsout),
            cbind(metric = rownames(nnet_out$confstatsout), nnet_out$confstatsout)
        )))
    colnames(modelstats_outtable) <- c("metric", paste0(models_ran, "_out"))
    ## Write out the modeling stats
    # write.table(modelstats_outtable, paste0(modeling_outfilepath, "modeling/model_stats_outtable.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    
    
    ## So this all works - the major issue is that the results from just using the test set actually looks ok... but just very blocky and off cause the N is so small. So we want to use the FULL data set, and to do that, we need to give the below data (predicted values and labels) for the full data set. Also doable, but we need to extract the appropriate predictor values per model.
    ## evalm aggregates using mean, by the rowindex of the BEST tuned model using the param listed in the model call/
    # 1 - grab best tune, 2 - grab rows from model WITH best tune, 3 - aggregate usingmean across Rowindexes for the subsetted best tune table, 4 - this data can then be input into roc
    # mm1 <- aggregate(model_rf[[5]][model_rf[[5]][,"mtry"] == 2,c(3,4)], list(model_rf[[5]][model_rf[[5]][,"mtry"] == 2,5]), mean)
    
    # rocdata <- cbind.data.frame(testing[,outcomelabel,drop=FALSE],
    #                             glm_out$predicted_values_and_labels,
    #                             xgbTree_out$predicted_values_and_labels,
    #                             svmRadial_out$predicted_values_and_labels,
    #                             rf_out$predicted_values_and_labels,
    #                             glmnet_out$predicted_values_and_labels)
    
    rocdata <- do.call(cbind.data.frame, Filter(Negate(is.null), list(
        testing[,outcomelabel,drop=FALSE],
        glm_out$predicted_values_and_labels,
        xgbTree_out$predicted_values_and_labels,
        svmRadial_out$predicted_values_and_labels,
        rf_out$predicted_values_and_labels,
        glmnet_out$predicted_values_and_labels,
        nnet_out$predicted_values_and_labels
    )))
    
    # Combine the data we want to use for plotting
    # rocdata <- cbind.data.frame(testing[,outcomelabel,drop=FALSE],
    #           predict_glm_values = predict_glm_values[,"sle"], predict_glm_labels = predict_glm_labels,
    #           predict_gbm_values = predict_gbm_values[,"sle"], predict_gbm_labels = predict_gbm_labels,
    #           predict_svmradial_values = predict_svmradial_values[,"sle"], predict_svmradial_labels = predict_svmradial_labels,
    #           predict_rf_values = predict_rf_values[,"sle"], predict_rf_labels = predict_rf_labels)
    
    # Combine the data we want to use for plotting
    # rocdata <- cbind.data.frame(testing[,outcomelabel,drop=FALSE], 
    #           predict_glm_values = predict_glm_values[,"sle"], predict_glm_labels = predict_glm_labels,
    #           predict_gbm_values = predict_gbm_values[,"sle"], predict_gbm_labels = predict_gbm_labels,
    #           predict_svmradial_values = predict_svmradial_values[,"sle"], predict_svmradial_labels = predict_svmradial_labels,
    #           predict_rf_values = predict_rf_values[,"sle"], predict_rf_labels = predict_rf_labels)
    
    
    
    ## Ok, so i need the true labels, and then the PREDICTOR value for SLE for each model, then if those are in one table, then we can run the below code and it should work.....
    # BUT we need the best tuend (or standard tuned) model values for each model which I think we can get using MLeval..
    
    # roc.list <- roc(diagnosis ~ predict_glm_values + predict_xgbTree_values + predict_svmRadial_values + predict_rf_values + predict_glmnet_values, 
    #                 data = rocdata, levels = rev(outcomelevels))
    # roc.list <- roc(as.formula(paste0(outcomelabel, "~ predict_glm_values + predict_xgbTree_values + predict_svmRadial_values + predict_rf_values + predict_glmnet_values")), data = rocdata, levels = rev(outcomelevels))
    # roc.list <- roc(as.formula(paste0(outcomelabel, "~ predict_svmRadial_values + predict_rf_values + predict_glmnet_values")), data = rocdata, levels = rev(outcomelevels))
    roc.list <- roc(as.formula(paste0(outcomelabel, "~ ", paste(colnames(rocdata)[grepl("_values", colnames(rocdata))], collapse = " + "))), data = rocdata, levels = rev(outcomelevels))
    if (class(roc.list) == "list") {
        ci_auc.list <- lapply(roc.list, ci.auc)
        auc_with_ci <- unlist(lapply(ci_auc.list, function(x) paste0("auc: ", round(x[[2]],2), " (", round(x[[1]],2), "-", round(x[[3]],2), ")")))
    } else {  ## if there is only one value
        ci_auc.list <- ci.auc(roc.list)
        # auc_with_ci <- unlist(lapply(ci_auc.list, function(x) paste0("auc: ", round(x[[2]],2), " (", round(x[[1]],2), "-", round(x[[3]],2), ")")))
        auc_with_ci <- paste0("auc: ", round(ci_auc.list[[2]],2), " (", round(ci_auc.list[[1]],2), "-", round(ci_auc.list[[3]],2), ")")
        names(auc_with_ci) <- colnames(rocdata)[2]
    }
    
    # auc_with_ci <- unlist(lapply(ci_auc.list, function(x) paste0("auc: ", round(x[[2]],2), " (", round(x[[1]],2), "-", round(x[[3]],2), ")")))
    auc_data_labels <- paste0(names(auc_with_ci), "\n", auc_with_ci)
    # lapply(roc.list, function(x) x[["auc"]])
    
    rocplot <- ggroc(roc.list) + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey")
    rocplot <- rocplot + scale_color_discrete(labels=auc_data_labels)
    rocplot <- rocplot + theme_minimal() + coord_equal()
    # pout
    
    # pdf(paste0(modeling_outfilepath, "modeling/ROCcurves.pdf"), 10, 10, useDingbats = FALSE)
    # print(rocplot)
    # junk <- dev.off()
    
    # Adds CI to ROC plots
    # ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))
    # dat.ci.list <- lapply(ci.list, function(ciobj) data.frame(x = as.numeric(rownames(ciobj)),
    #              lower = ciobj[, 1], upper = ciobj[, 3]))
    # for(i in seq_len(length(ci.list))) {
    #   p <- p + geom_ribbon(data = dat.ci.list[[i]], aes(x = x, ymin = lower, ymax = upper),
    #     fill = i + 1, alpha = 0.2, inherit.aes = F)
    # }

    return(list(modelstats_outtable=rbind(modelstats_outtable, c("AUC_w_CI", auc_with_ci)),
                rocplot=rocplot))
    
    
}





# --------------------------------- regression ML pipeline --------------------------------- 
regression_ml_analysis <- function(featuretable, outcometable, outcomelevels, seedparam=11111, train_partition_percent = 0.6,
                                              # subsampleparam = NULL, 
                                   models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")) {
    ## ML for SLE using Genes
    # packagelist = c("mlbench", "caret", "pROC")
    # junk <- lapply(packagelist, function(xxx) suppressMessages(
    #   require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
    # dir.create(paste0(modeling_outfilepath, "modeling/"), showWarnings = FALSE, recursive = TRUE)
    
    # Create modeling input data
    outcomelabel <- colnames(outcometable)
    featurelabels <- colnames(featuretable)
    modeling_intable <- na.omit(merge(outcometable, featuretable, by = "row.names"))
    rownames(modeling_intable) <- modeling_intable[,"Row.names"]
    modeling_intable <- modeling_intable[,!grepl("Row.names", colnames(modeling_intable))]
    # modeling_intable[,outcomelabel] <- factor(modeling_intable[,outcomelabel], levels = outcomelevels)
    
    # Ok - we need to adjust variable names where applicable (like with 3v and 1v - doesnt like numbers first)
    # level_key <- setNames(make.names(make.names(outcomelevels)), outcomelevels)
    # modeling_intable[,outcomelabel] <- recode(modeling_intable[,outcomelabel], !!!level_key)
    # outcomelevels <- make.names(outcomelevels)
    
    
    # create data partition and prepare training scheme
    #1415 0.92
    set.seed(seedparam)
    trainindex <- createDataPartition(modeling_intable[,outcomelabel], p = train_partition_percent, list = FALSE)
    training <- modeling_intable[trainindex,]
    testing <- modeling_intable[-trainindex,]

    
    # trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, savePredictions = TRUE)
    # trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5)
    trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5)

    
    regression_modeltraining <- function(training, testing, outcomelabel, featurelabels, trcontrol, modellabel, seedparam) {
        set.seed(seedparam)
        print(paste0("Training ", modellabel, " model"))
        modfit <- train(as.formula(paste0(outcomelabel, "~ .")), data = training, method = modellabel, 
                        trControl=trcontrol, preProc = c("center", "scale"))
        if (modellabel == "xgbTree") {
            modfit <- train(as.formula(paste0(outcomelabel, "~ .")), data = training, method = modellabel, 
                            # trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC", verbose = FALSE, verbosity = 0)
                            trControl=trcontrol, preProc = c("center", "scale"), verbose = FALSE, verbosity = 0)
        }
        # predict_labels <- predict(modfit, testing[,featurelabels])
        predict_values <- predict(modfit, testing)
        # plot(testing[,outcomelabel], predict_values)
        # predict_values <- predict(modfit, testing, type = "prob")[,1,drop=FALSE]
        # confstatsout <- writeout_confmat(confusionMatrix(predict_labels, testing[,outcomelabel]))
        predicted_values_and_labels <- cbind(predict_values, testing[,outcomelabel])
        colnames(predicted_values_and_labels) <- c(paste0("predict_", modellabel, "_values"), paste0("actual_", modellabel, "_values"))
        return(list(predicted_values_and_labels = predicted_values_and_labels))
    }
    # train some models
    
    # Im gonna do them separately instead of a for loop for now - probably not as neat, but easier to track for now
    models_ran <- c()
    if ("glm" %in% models_to_run) {
        glm_out <- regression_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                              modellabel = "glm", seedparam = seedparam)
        models_ran <- c(models_ran, "glm")
    } else { glm_out <- NULL }
    if ("xgbTree" %in% models_to_run) {
        xgbTree_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                  modellabel = "xgbTree", seedparam = seedparam)
        models_ran <- c(models_ran, "xgbTree")
    } else { xgbTree_out <- NULL }
    if ("svmRadial" %in% models_to_run) {
        svmRadial_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                    modellabel = "svmRadial", seedparam = seedparam)
        models_ran <- c(models_ran, "svmRadial")
    } else { svmRadial_out <- NULL }
    if ("rf" %in% models_to_run) {
        rf_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                             modellabel = "rf", seedparam = seedparam)
        models_ran <- c(models_ran, "rf")
    } else { rf_out <- NULL }
    if ("glmnet" %in% models_to_run) {
        glmnet_out <- cohort_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                 modellabel = "glmnet", seedparam = seedparam) ## lasso
        models_ran <- c(models_ran, "glmnet")
    } else { glmnet_out <- NULL }
    
    
    ## Attach to inmetatable by merge (for safety) then reassign and remove rownames
    modelstats_outtable = Reduce(function(dtf1, dtf2)
        merge(dtf1, dtf2, by = "metric", all.x = TRUE, sort = FALSE),
        Filter(Negate(is.null), list(
            cbind(metric = rownames(glm_out$confstatsout), glm_out$confstatsout),
            cbind(metric = rownames(xgbTree_out$confstatsout), xgbTree_out$confstatsout),
            cbind(metric = rownames(svmRadial_out$confstatsout), svmRadial_out$confstatsout), 
            cbind(metric = rownames(rf_out$confstatsout), rf_out$confstatsout), 
            cbind(metric = rownames(glmnet_out$confstatsout), glmnet_out$confstatsout) 
        )))
    colnames(modelstats_outtable) <- c("metric", paste0(models_ran, "_out"))
    ## Write out the modeling stats
    # write.table(modelstats_outtable, paste0(modeling_outfilepath, "modeling/model_stats_outtable.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    
    
    ## So this all works - the major issue is that the results from just using the test set actually looks ok... but just very blocky and off cause the N is so small. So we want to use the FULL data set, and to do that, we need to give the below data (predicted values and labels) for the full data set. Also doable, but we need to extract the appropriate predictor values per model.
    ## evalm aggregates using mean, by the rowindex of the BEST tuned model using the param listed in the model call/
    # 1 - grab best tune, 2 - grab rows from model WITH best tune, 3 - aggregate usingmean across Rowindexes for the subsetted best tune table, 4 - this data can then be input into roc
    # mm1 <- aggregate(model_rf[[5]][model_rf[[5]][,"mtry"] == 2,c(3,4)], list(model_rf[[5]][model_rf[[5]][,"mtry"] == 2,5]), mean)
    
    # rocdata <- cbind.data.frame(testing[,outcomelabel,drop=FALSE],
    #                             glm_out$predicted_values_and_labels,
    #                             xgbTree_out$predicted_values_and_labels,
    #                             svmRadial_out$predicted_values_and_labels,
    #                             rf_out$predicted_values_and_labels,
    #                             glmnet_out$predicted_values_and_labels)
    
    rocdata <- do.call(cbind.data.frame, Filter(Negate(is.null), list(
        testing[,outcomelabel,drop=FALSE],
        glm_out$predicted_values_and_labels,
        xgbTree_out$predicted_values_and_labels,
        svmRadial_out$predicted_values_and_labels,
        rf_out$predicted_values_and_labels,
        glmnet_out$predicted_values_and_labels
    )))
    
    # Combine the data we want to use for plotting
    # rocdata <- cbind.data.frame(testing[,outcomelabel,drop=FALSE],
    #           predict_glm_values = predict_glm_values[,"sle"], predict_glm_labels = predict_glm_labels,
    #           predict_gbm_values = predict_gbm_values[,"sle"], predict_gbm_labels = predict_gbm_labels,
    #           predict_svmradial_values = predict_svmradial_values[,"sle"], predict_svmradial_labels = predict_svmradial_labels,
    #           predict_rf_values = predict_rf_values[,"sle"], predict_rf_labels = predict_rf_labels)
    
    # Combine the data we want to use for plotting
    # rocdata <- cbind.data.frame(testing[,outcomelabel,drop=FALSE], 
    #           predict_glm_values = predict_glm_values[,"sle"], predict_glm_labels = predict_glm_labels,
    #           predict_gbm_values = predict_gbm_values[,"sle"], predict_gbm_labels = predict_gbm_labels,
    #           predict_svmradial_values = predict_svmradial_values[,"sle"], predict_svmradial_labels = predict_svmradial_labels,
    #           predict_rf_values = predict_rf_values[,"sle"], predict_rf_labels = predict_rf_labels)
    
    
    
    ## Ok, so i need the true labels, and then the PREDICTOR value for SLE for each model, then if those are in one table, then we can run the below code and it should work.....
    # BUT we need the best tuend (or standard tuned) model values for each model which I think we can get using MLeval..
    
    # roc.list <- roc(diagnosis ~ predict_glm_values + predict_xgbTree_values + predict_svmRadial_values + predict_rf_values + predict_glmnet_values, 
    #                 data = rocdata, levels = rev(outcomelevels))
    # roc.list <- roc(as.formula(paste0(outcomelabel, "~ predict_glm_values + predict_xgbTree_values + predict_svmRadial_values + predict_rf_values + predict_glmnet_values")), data = rocdata, levels = rev(outcomelevels))
    # roc.list <- roc(as.formula(paste0(outcomelabel, "~ predict_svmRadial_values + predict_rf_values + predict_glmnet_values")), data = rocdata, levels = rev(outcomelevels))
    roc.list <- roc(as.formula(paste0(outcomelabel, "~ ", paste(colnames(rocdata)[grepl("_values", colnames(rocdata))], collapse = " + "))), data = rocdata, levels = rev(outcomelevels))
    ci_auc.list <- lapply(roc.list, ci.auc)
    auc_with_ci <- unlist(lapply(ci_auc.list, function(x) paste0("auc: ", round(x[[2]],2), " (", round(x[[1]],2), "-", round(x[[3]],2), ")")))
    auc_data_labels <- paste0(names(auc_with_ci), "\n", auc_with_ci)
    # lapply(roc.list, function(x) x[["auc"]])
    
    rocplot <- ggroc(roc.list) + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey")
    rocplot <- rocplot + scale_color_discrete(labels=auc_data_labels)
    rocplot <- rocplot + theme_minimal() + coord_equal()
    # pout
    
    # pdf(paste0(modeling_outfilepath, "modeling/ROCcurves.pdf"), 10, 10, useDingbats = FALSE)
    # print(rocplot)
    # junk <- dev.off()
    
    # Adds CI to ROC plots
    # ci.list <- lapply(roc.list, ci.se, specificities = seq(0, 1, l = 25))
    # dat.ci.list <- lapply(ci.list, function(ciobj) data.frame(x = as.numeric(rownames(ciobj)),
    #              lower = ciobj[, 1], upper = ciobj[, 3]))
    # for(i in seq_len(length(ci.list))) {
    #   p <- p + geom_ribbon(data = dat.ci.list[[i]], aes(x = x, ymin = lower, ymax = upper),
    #     fill = i + 1, alpha = 0.2, inherit.aes = F)
    # }
    
    return(list(modelstats_outtable=rbind(modelstats_outtable, c("AUC_w_CI", auc_with_ci)),
                rocplot=rocplot))
    
    
}


# --------------------------------- regression TO clasifiation ML pipeline --------------------------------- 

regression_to_classification_ml_analysis <- function(featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
                                                      seedparam=11111, subsampleparam = NULL, train_partition_percent = 0.6,
                                                      models_to_run = c("glm", "svmRadial", "rf", "glmnet")) {
    # Create modeling input data
    outcomelabel_time <- colnames(outcometable)[grepl("T_", colnames(outcometable))]
    outcomelabel_class <- colnames(outcometable)[grepl("C_", colnames(outcometable))]
    featurelabels <- colnames(featuretable)
    modeling_intable <- na.omit(merge(outcometable, featuretable, by = "row.names"))
    rownames(modeling_intable) <- modeling_intable[,"Row.names"]
    modeling_intable <- modeling_intable[,!grepl("Row.names", colnames(modeling_intable))]
    
    # create data partition and prepare training scheme
    set.seed(seedparam)
    trainindex <- createDataPartition(modeling_intable[,outcomelabel_class], p = train_partition_percent, list = FALSE)
    reg_training <- modeling_intable[trainindex, c(outcomelabel_time, featurelabels)]
    reg_testing <- modeling_intable[-trainindex, c(outcomelabel_time, featurelabels)]
    class_training <- modeling_intable[trainindex, c(outcomelabel_class, featurelabels)]
    class_testing <- modeling_intable[-trainindex, c(outcomelabel_class, featurelabels)]

    # Now first - run our TRAIN DATA through a REGRESSION analysis
    print("training the time to event prediction - regression")
    trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, savePredictions = )
    regression_modeltraining <- function(training, testing, outcomelabel, featurelabels, trcontrol, modellabel, seedparam) {
        set.seed(seedparam)
        print(paste0("Training ", modellabel, " model"))
        modfit <- train(as.formula(paste0(outcomelabel, "~ .")), data = training, method = modellabel, 
                        trControl=trcontrol, preProc = c("center", "scale"))
        if (modellabel == "xgbTree") {
            modfit <- train(as.formula(paste0(outcomelabel, "~ .")), data = training, method = modellabel, 
                            trControl=trcontrol, preProc = c("center", "scale"), verbose = FALSE, verbosity = 0)
        }
        # Predict on our testing
        predict_values <- predict(modfit, testing)
        predicted_values_and_labels <- cbind(predict_values, testing[,outcomelabel,drop=FALSE], sampleID = rownames(testing))
        colnames(predicted_values_and_labels) <- c(paste0("predict_", modellabel, "_values"), paste0("actual_", modellabel, "_values"), "sampleID")
        
        # predict_values <- predict(modfit, training)
        # predicted_values_and_labels <- cbind(predict_values, training[,outcomelabel,drop=FALSE], sampleID = rownames(training))

        return(list(predicted_values_and_labels = predicted_values_and_labels, modelout = modfit))
    }
    # train some models
    
    # for each model - grab the predicted outputs
    test_pred_values_out <- train_pred_values_out <- list()
    for (model_sel in models_to_run) {
        reg_training_out <- regression_modeltraining(training=reg_training, testing=reg_testing,
                                                     outcomelabel=outcomelabel_time, featurelabels=featurelabels, trcontrol,
                                                     modellabel=model_sel, seedparam=seedparam)
        # grab the predicted values for test:
        predition_table_test <- reg_training_out[["predicted_values_and_labels"]]
        # colnames(predition_table)[2] <- "actual_values"
        test_pred_values_out[[model_sel]] <- predition_table_test
        
        # We dont actually want the TEST times out - we want the TRAIN times to TRAIN our class model
        ## So use the trained model to grab our predicted values for the TRAIN set
        modfit <- reg_training_out$modelout
        predition_table_train <- cbind(sampleID = rownames(reg_training[,outcomelabel_time,drop=FALSE]), 
                                 data.frame(predict(modfit, reg_training)), reg_training[,outcomelabel_time,drop=FALSE])
        colnames(predition_table_train) <- c("sampleID", paste0("predict_", model_sel, "_values"), "actual_values")
        train_pred_values_out[[model_sel]] <- predition_table_train
        
        
    }
    test_prediction_table_full <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sampleID", all = TRUE, sort = FALSE), test_pred_values_out)
    rownames(test_prediction_table_full) <- test_prediction_table_full[,"sampleID"]
    test_prediction_table_full <- test_prediction_table_full[,grepl("predict_", colnames(test_prediction_table_full))]
    train_prediction_table_full <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sampleID", all = TRUE, sort = FALSE), train_pred_values_out)
    rownames(train_prediction_table_full) <- train_prediction_table_full[,"sampleID"]
    train_prediction_table_full <- train_prediction_table_full[,grepl("predict_", colnames(train_prediction_table_full))]
    
    # Now need to train our class model - but i dont know if I should use the ORIGINAL time to event data to TRAIN, or our NEW data
    ## Regardless - for testing - we use our OUTPUTTED DATA from our reg train/test (train_prediction_table_full)
    
    # And regardless - we still have our trains/split - so we need to keep the train split
    # for each model - grab the predicted outputs
    class_stat_out <- list()
    print("training the event prediction - classification")
    for (model_sel in models_to_run) {
        # Grab one of our reg model outputs for INPUT into this table
        featuretable <- train_prediction_table_full[,grepl(paste0("_", model_sel, "_"), colnames(train_prediction_table_full)),drop=FALSE]
        outcometable <- class_training[,1,drop=FALSE]
        presplit_training <- cbind(outcometable, featuretable)
        presplit_testing <- cbind(class_testing[,1,drop=FALSE], 
                                  test_prediction_table_full[,grepl(paste0("_", model_sel, "_"), colnames(test_prediction_table_full)),drop=FALSE])
        # Then run our modeling with these genes and labels
        class_models_to_run <- models_to_run
        if ("glmnet" %in% models_to_run) {
            class_models_to_run <- class_models_to_run[!class_models_to_run %in% "glmnet"]
            print("cant run glmnet with one feature (time) so we are excluding")
        }
        presplitdata <- list(training = presplit_training, testing = presplit_testing)
        out1 <- cohort_classification_ml_analysis(featuretable = "dummy", outcometable = "dummy", 
                                                  outcomelevels = outcomelevels,seedparam=seedparam, presplitdata = presplitdata,
                                                  subsampleparam = subsampleparam, models_to_run = class_models_to_run)
        # write out results
        modelstatout <- out1[[1]]
        colnames(modelstatout)[2:ncol(modelstatout)] <- c(paste0(model_sel, "_reg__", colnames(modelstatout)[2:ncol(modelstatout)]))
        class_stat_out[[model_sel]] <- modelstatout
    }
    class_stat_table_out <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "metric", all = TRUE, sort = FALSE), class_stat_out)
    return(class_stat_table_out)
}




# --------------------------------- sample scoring for using singscore and genesets --------------------------------- 

singscore_sample_scoring <- function(GOIup, GOIdown, counttable) {
    ## Trying singscore
    library(singscore)
    
    ## first rank genes with the rankgene function
    rankData <- rankGenes(counttable)
    # Then calculate the score
    scoredf <- simpleScore(rankData, upSet = GOIup, downSet = GOIdown)
    
    return(scoredf)
    
}

# --------------------------------- sample socring for WGCNA eigengenes --------------------------------- 
WGNCA_sample_scoring <- function(GOIup, GOIdown, counttable, complabel_sel) {
    library("WGCNA")
    
    # transpose to make WGCNA format
    WGCNA_counttable <- t(counttable)
    
    # Run for up, down, and both (i guess?)
    # genestocolors <- data.frame(moduleLabels = ifelse(colnames(WGCNA_counttable) %in% GOIup, 1, 0), 
    #                             moduleColors = ifelse(colnames(WGCNA_counttable) %in% GOIup, complabel_sel, "grey"), 
    #                             row.names = colnames(WGCNA_counttable), stringsAsFactors = FALSE)
    # moduleColors <- genestocolors[,"moduleColors"]
    
    # Run for all GOI
    moduleColors_all <- ifelse(colnames(WGCNA_counttable) %in% c(GOIup, GOIdown), "WGCNAscore_all", "grey")
    eigengenes_all <- orderMEs(moduleEigengenes(WGCNA_counttable, moduleColors_all)$eigengenes)
    
    # Run for up GOI
    moduleColors_up <- ifelse(colnames(WGCNA_counttable) %in% GOIup, "WGCNAscore_up", "grey")
    eigengenes_up <- orderMEs(moduleEigengenes(WGCNA_counttable, moduleColors_up)$eigengenes)
    
    # Run for down GOI
    moduleColors_down <- ifelse(colnames(WGCNA_counttable) %in% GOIdown, "WGCNAscore_down", "grey")
    eigengenes_down <- orderMEs(moduleEigengenes(WGCNA_counttable, moduleColors_down)$eigengenes)
    
    # Combine all to write out
    eigengene_fulltable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "SampleID", all = TRUE, sort = FALSE),
                                  list(cbind.data.frame(SampleID = rownames(eigengenes_all), eigengenes_all),
                                       cbind.data.frame(SampleID = rownames(eigengenes_up), eigengenes_up),
                                       cbind.data.frame(SampleID = rownames(eigengenes_down), eigengenes_down)))
    rownames(eigengene_fulltable) <- eigengene_fulltable[,1]
    colnames(eigengene_fulltable) <- c("SampleID", "WGCNAscore_all", "MEgrey_all", "WGCNAscore_up", "MEgrey_up", "WGCNAscore_down", "MEgrey_down")
    
    return(eigengene_fulltable[,!grepl("SampleID", colnames(eigengene_fulltable))])
    
}

# --------------------------------- signature scoring for samples using both methods wrap up function --------------------------------- 
sample_signature_analysis <- function(intable_deseq, dge_comp_label, pval_test_sel, pval_cutoff_sel, log2fc_cutoff_sel,
                                      comp_metatable, counttable_sel, signaturescore_outfilepath, run_singscore = TRUE, run_WGCNAscore = TRUE,
                                      scoring_genesetsize = c("all", 1000, 500, 100, 50)) {
    # GOI UP
    # Create a ranking metric for these genes, perhaps deciles of 3 metrics: log2fc, padj, and baseMean - then take sum of deciles?
    GOIup_table <- intable_deseq[intable_deseq[,pval_test_sel] < pval_cutoff_sel & !is.na(intable_deseq[,pval_test_sel]) & 
                                     intable_deseq[,"log2FoldChange"] > log2fc_cutoff_sel, ]
    GOIup_table[,"pstat"] <- -log10(GOIup_table[,pval_test_sel])
    GOIup_table[,"abslog2FoldChange"] <- abs(GOIup_table[,"log2FoldChange"])
    
    scoringmetrics <- c("baseMean", "abslog2FoldChange", "pstat")[c("baseMean", "abslog2FoldChange", "pstat") %in% colnames(GOIup_table)]
    
    out1 <- apply(GOIup_table[,scoringmetrics], 2, function(x) 
        continuous_to_named_quantile(data.frame(x), number_of_quantiles = 10))
    # From this - I need to grab the quantiles from each (for reference), combine the quantile tables and output those
    quantile_grab <- do.call(cbind, lapply(out1, function(x) as.character(x[["quantiles"]])))
    quantiletable_grab <- do.call(cbind, lapply(out1, function(x) {
        quant_df <- x[["quantile_dataframe"]]
        quant_df[,1] <- as.numeric(gsub("Q|_Low|_High", "", as.character(quant_df[,1])))
        quant_df
    }))
    colnames(quantiletable_grab) <- scoringmetrics
    quantiletable_grab[,"quant_sum"] <- rowSums(quantiletable_grab[,scoringmetrics])
    
    GOIup_table_wstats <- cbind(GOIup_table, quantiletable_grab)
    GOIup_table_wstats <- GOIup_table_wstats[order(GOIup_table_wstats[,"quant_sum"], GOIup_table_wstats[,"abslog2FoldChange"], decreasing = TRUE),]
    GOIup <- rownames(GOIup_table_wstats)
    
    write.table(quantile_grab, paste0(signaturescore_outfilepath, dge_comp_label, "_GOIup_quantiles.csv"),
                sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(data.frame(GOIup_table_wstats), paste0(signaturescore_outfilepath, dge_comp_label, "_GOIup_table_wstats.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    write.table(data.frame(GOIup), paste0(signaturescore_outfilepath, dge_comp_label, "_GOIup.csv"),
                sep = ",", col.names = FALSE, row.names = FALSE)
    
    # GOI down
    GOIdown_table <- intable_deseq[intable_deseq[,pval_test_sel] < pval_cutoff_sel & !is.na(intable_deseq[,pval_test_sel]) & 
                                     intable_deseq[,"log2FoldChange"] < -log2fc_cutoff_sel, ]
    GOIdown_table[,"pstat"] <- -log10(GOIdown_table[,pval_test_sel])
    GOIdown_table[,"abslog2FoldChange"] <- abs(GOIdown_table[,"log2FoldChange"])
    out1 <- apply(GOIdown_table[,scoringmetrics], 2, function(x) 
        continuous_to_named_quantile(data.frame(x), number_of_quantiles = 10))
    # From this - I need to grab the quantiles from each (for reference), combine the quantile tables and output those
    quantile_grab <- do.call(cbind, lapply(out1, function(x) as.character(x[["quantiles"]])))
    quantiletable_grab <- do.call(cbind, lapply(out1, function(x) {
        quant_df <- x[["quantile_dataframe"]]
        quant_df[,1] <- as.numeric(gsub("Q|_Low|_High", "", as.character(quant_df[,1])))
        quant_df
    }))
    colnames(quantiletable_grab) <- scoringmetrics
    quantiletable_grab[,"quant_sum"] <- rowSums(quantiletable_grab[,scoringmetrics])
    
    GOIdown_table_wstats <- cbind(GOIdown_table, quantiletable_grab)
    GOIdown_table_wstats <- GOIdown_table_wstats[order(GOIdown_table_wstats[,"quant_sum"], GOIdown_table_wstats[,"abslog2FoldChange"], decreasing = TRUE),]
    GOIdown <- rownames(GOIdown_table_wstats)
    
    write.table(quantile_grab, paste0(signaturescore_outfilepath, dge_comp_label, "_GOIdown_quantiles.csv"),
                sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(data.frame(GOIdown_table_wstats), paste0(signaturescore_outfilepath, dge_comp_label, "_GOIdown_table_wstats.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    write.table(data.frame(GOIdown), paste0(signaturescore_outfilepath, dge_comp_label, "_GOIdown.csv"),
                sep = ",", col.names = FALSE, row.names = FALSE)
    
    for (genesetsize in scoring_genesetsize) {
        if (genesetsize == "all") {
            select_GOIup <- GOIup
            select_GOIdown <- GOIdown
            UPgenesetsize <- DOWNgenesetsize <- "all"
        } else {  ## Is a number
            UPgenesetsize <- min(nrow(GOIup_table_wstats), as.numeric(genesetsize))
            DOWNgenesetsize <- min(nrow(GOIdown_table_wstats), as.numeric(genesetsize))
            
            select_GOIup <- rownames(GOIup_table_wstats)[1:UPgenesetsize]
            select_GOIdown <- rownames(GOIdown_table_wstats)[1:DOWNgenesetsize]
        }
        
        # sing score first
        if (run_singscore) {
            subanalysis_outfilepath <- paste0(signaturescore_outfilepath, "singscore/", UPgenesetsize, "UP_", DOWNgenesetsize, "DOWN_genes/")
            dir.create(subanalysis_outfilepath, showWarnings = FALSE, recursive = TRUE)
            
            # generate the score
            singscore_table <- singscore_sample_scoring(GOIup=select_GOIup, GOIdown=select_GOIdown, counttable=counttable_sel)
            # ggpairs(singscore_table[,c("TotalScore", "UpScore", "DownScore")])
            
            # Now what should we output...
            # (1) the score table itself, (2) the gene sets, (3) I really think a boxplot for all of our samples using the metric is good for each one
            singscore_outtable <- merge(comp_metatable, singscore_table, by = "row.names")
            write.table(singscore_outtable, paste0(subanalysis_outfilepath, dge_comp_label, "_singscore_outtable.csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            write.table(select_GOIup, paste0(subanalysis_outfilepath, dge_comp_label, "_GOIup_", UPgenesetsize, ".csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            write.table(select_GOIdown, paste0(subanalysis_outfilepath, dge_comp_label, "_GOIdown_", DOWNgenesetsize, ".csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            
            
            ## Plot a box plot just to check if it works
            for (compmetacolnum in seq_len(ncol(comp_metatable))) {
                compmetacol_sel <- colnames(comp_metatable[,compmetacolnum,drop=FALSE])
                pdf(paste0(subanalysis_outfilepath, compmetacol_sel, "_singscore_bpplot.pdf"))
                for (scoretype in c("TotalScore", "TotalDispersion", "UpScore", "UpDispersion", "DownScore", "DownDispersion")) {
                    bptable <- singscore_outtable[,c("Row.names", compmetacol_sel, scoretype)]
                    bptable[,1] <- compmetacol_sel
                    bpout1 <- boxplot_plotter(boxplottable = na.omit(bptable), xsplit = "feature", plotstats = "intra", 
                                              labsparam = list(title = paste0(compmetacol_sel, "__singscore_", scoretype), x = compmetacol_sel, y = scoretype,
                                                               catorder = sort(unique(bptable[,compmetacol_sel])), featorder = compmetacol_sel))
                    
                    print(bpout1)
                }
                junk <- dev.off()
            }
        }
        
        
        # then WGCNA
        if (run_WGCNAscore) {
            subanalysis_outfilepath <- paste0(signaturescore_outfilepath, "WGCNAscore/", UPgenesetsize, "UP_", DOWNgenesetsize, "DOWN_genes/")
            dir.create(subanalysis_outfilepath, showWarnings = FALSE, recursive = TRUE)
            # generate the score
            WGCNAscore_table <- WGNCA_sample_scoring(GOIup=select_GOIup, GOIdown=select_GOIdown, counttable=counttable_sel, complabel_sel=complabel_sel)
            # ggpairs(singscore_table[,c("TotalScore", "UpScore", "DownScore")])
            
            # Now what should we output...
            # (1) the score table itself, (2) the gene sets, (3) I really think a boxplot for all of our samples using the metric is good for each one
            WGCNAscore_outtable <- merge(comp_metatable, WGCNAscore_table, by = "row.names")
            write.table(WGCNAscore_outtable, paste0(subanalysis_outfilepath, dge_comp_label, "_WGCNAscore_outtable.csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            write.table(select_GOIup, paste0(subanalysis_outfilepath, dge_comp_label, "_GOIup_", UPgenesetsize, ".csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            write.table(select_GOIdown, paste0(subanalysis_outfilepath, dge_comp_label, "_GOIdown_", DOWNgenesetsize, ".csv"),
                        sep = ",", col.names = TRUE, row.names = FALSE)
            
            
            ## Plot a box plot just to check if it works
            for (compmetacolnum in seq_len(ncol(comp_metatable))) {
                compmetacol_sel <- colnames(comp_metatable[,compmetacolnum,drop=FALSE])
                bptable <- WGCNAscore_outtable[,c("Row.names", compmetacol_sel, "WGCNAscore_all")]
                bptable[,1] <- compmetacol_sel
                bpout1 <- boxplot_plotter(boxplottable = na.omit(bptable), xsplit = "feature", plotstats = "intra", 
                                          labsparam = list(title = paste0(compmetacol_sel, "__WGCNAscore"), x = compmetacol_sel, y = "WGCNAscore_all",
                                                           catorder = sort(unique(bptable[,compmetacol_sel])), featorder = compmetacol_sel))
                pdf(paste0(subanalysis_outfilepath, compmetacol_sel, "_WGCNAscore_bpplot.pdf"))
                print(bpout1)
                junk <- dev.off()
            }
        }
    }
    
}

# --------------------------------- filter a time to event table into binary by time filters --------------------------------- 
# Filters a time to event table by the specified time filters to return events within a certain time
filter_time_to_event_table <- function(time_to_event_table, time_filters) {
    for (timevar in time_filters) {
        EOI <- gsub("T_", "", colnames(time_to_event_table)[1])
        time_to_event_table[,paste0("filtered", timevar, "_", EOI)] <- 
            ifelse(time_to_event_table[, paste0("T_", EOI)] >= timevar & time_to_event_table[, paste0("C_", EOI)] == 1, 2, # filter out too late event
            ifelse(time_to_event_table[, paste0("T_", EOI)] < timevar & time_to_event_table[, paste0("C_", EOI)] == 0, 2, # filter out too soon non-event
                   time_to_event_table[, paste0("C_", EOI)]))
    }
    return(time_to_event_table)
}



# --------------------------------- quantiles a Nx1 dataframe --------------------------------- 
# Turns a Nx1 data frame into a quantiled dataframe
continuous_to_named_quantile <- function(continuous_dataframe, number_of_quantiles) {
    cut_dataframe <- continuous_dataframe
    cut_dataframe[,1] <- cut2(t(continuous_dataframe), g=number_of_quantiles)
    original_ref_dataframe <- cbind(continuous_dataframe, cut2(t(continuous_dataframe), g=number_of_quantiles))
    quantile_ranges <- levels(cut_dataframe[,1])
    for (quantile_num in seq_len(length(quantile_ranges))) {
        quantile_sel <- quantile_ranges[quantile_num]
        quantile_label <- paste0("Q", quantile_num)
        if (quantile_num == 1) {quantile_label <- paste0(quantile_label, "_Low")}
        if (quantile_num == length(quantile_ranges)) {quantile_label <- paste0(quantile_label, "_High")}
        levels(cut_dataframe[,1])[match(levels(cut_dataframe[,1])[quantile_num], levels(cut_dataframe[,1]))] <- quantile_label
    }
    return(list(quantile_dataframe=cut_dataframe, original_ref_dataframe=original_ref_dataframe, quantiles = unique(sort(original_ref_dataframe[,2]))))
}


# --------------------------------- Deconvolution function --------------------------------- 

# Deconvolution analysis
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(readr)
library(parallel)
library(e1071)
library(preprocessCore)
set_cibersort_binary("/Users/tosh/Desktop/Ruggles_Lab/code/immunedeconv_scripts/cibersort_files/CIBERSORT_MGC_CUSTOM.R")
set_cibersort_mat("/Users/tosh/Desktop/Ruggles_Lab/code/immunedeconv_scripts/cibersort_files/LM22.txt")

# MCPcounter             EPIC        quanTIseq            xCell        CIBERSORT CIBERSORT (abs.)  TIMER 
# "mcp_counter"           "epic"      "quantiseq"          "xcell"      "cibersort"  "cibersort_abs"   "timer" 
# deconv_alg <- c("mcp_counter", "epic", "quantiseq", "xcell", "cibersort", "cibersort_abs")

deconvolution_analysis <- function(normcounttable, deconv_alg) {
    deconv_table_outlist <- list()
    for (alg_sel in deconv_alg) {
        deconv_res <- deconvolute(normcounttable, alg_sel)
        # Reformat the cibersort output
        rownames(deconv_res) <- paste0(alg_sel, "_", unlist(deconv_res[,"cell_type"]))
        deconv_outtable <- data.frame(t(deconv_res)[!rownames(t(deconv_res)) %in% "cell_type",])
        deconv_outtable[] <- apply(deconv_outtable, 2, as.numeric)
        
        deconv_table_outlist[[alg_sel]] <- deconv_outtable
    }
    return(deconv_table_outlist)
}




# Manual ROC curve for a 2-axis analysis
# Input the table with the 2 metric columns and the classifier
# Will also need the levels of the classifier AND the metics if they are ordinal
# Then it will basically plot the scatterplot (easy) and then TRACK along the diagonal (but which diagonal.......)
# On second thought - think itll have to be just confusion matrices with accuracies at different cut offs - then a grid search for the best values
# so along those lines, we need a grid of cut offs to use for the plot...... And those will be 

# predictor_table <- bioreptable_waddons[,c("CUSTOM_IMGDEGIS_NONEMILD", "CTNDV70", "")]









# --------------------------------- cohort multi class ml function ---------------------------------

cohort_multi_classification_ml_analysis <- function(featuretable, outcometable, outcomelevels, seedparam=11111,
                                                    presplitdata = NULL, train_partition_percent = 0.6, OHE_featuretable = TRUE,
                                                    subsampleparam = NULL, models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet", "nnet", "kknn")) {
    
    # One hot encoding of categorical features
    if (OHE_featuretable) {
        featuretable <- data.frame(predict(dummyVars(~., featuretable), featuretable))
    }
    
    # xx1 <- data.frame(predict(dummyVars(~., featuretable), featuretable))
    # xx2 <- xx1[,!colnames(xx1) %in% colnames(featuretable)]
    # t(apply(xx2, 2, table))
    
    # Create modeling input data
    if (is.null(presplitdata)) {
        outcomelabel <- colnames(outcometable)
        featurelabels <- colnames(featuretable)
        modeling_intable <- na.omit(merge(outcometable, featuretable, by = "row.names"))
        rownames(modeling_intable) <- modeling_intable[,"Row.names"]
        modeling_intable <- modeling_intable[,!grepl("Row.names", colnames(modeling_intable))]
        modeling_intable[,outcomelabel] <- factor(modeling_intable[,outcomelabel], levels = outcomelevels)
        
        # Ok - we need to adjust variable names where applicable (like with 3v and 1v - doesnt like numbers first)
        level_key <- setNames(make.names(make.names(outcomelevels)), outcomelevels)
        modeling_intable[,outcomelabel] <- recode(modeling_intable[,outcomelabel], !!!level_key)
        outcomelevels <- make.names(outcomelevels)
        
        
        # create data partition and prepare training scheme
        #1415 0.92
        set.seed(seedparam)
        trainindex <- createDataPartition(modeling_intable[,outcomelabel], p = train_partition_percent, list = FALSE)
        training <- modeling_intable[trainindex,]
        testing <- modeling_intable[-trainindex,]
    } else {
        training <- presplitdata[["training"]]
        testing <- presplitdata[["testing"]]
        outcomelabel <- colnames(training)[1]
        featurelabels <- colnames(training)[2:ncol(training)]
        
        # Still have to do the name catch/change - for training and testing separately
        level_key <- setNames(make.names(make.names(outcomelevels)), outcomelevels)
        training[,outcomelabel] <- recode(training[,outcomelabel], !!!level_key)
        testing[,outcomelabel] <- recode(testing[,outcomelabel], !!!level_key)
        outcomelevels <- make.names(outcomelevels)
        training[,outcomelabel] <- factor(training[,outcomelabel], levels = outcomelevels)
        testing[,outcomelabel] <- factor(testing[,outcomelabel], levels = outcomelevels)
    }
    # trainindex <- createDataPartition(modeling_intable[,outcomelabel], p = 1, list = FALSE)
    # training <- modeling_intable[trainindex,]
    # testing <- modeling_intable[trainindex,]
    if (is.null(subsampleparam)) {
        trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, 
                                  classProbs = TRUE, savePredictions = TRUE, summaryFunction = multiClassSummary)
    } else {
        trcontrol <- trainControl(method="repeatedcv", number=10, repeats=5, 
                                  classProbs = TRUE, savePredictions = TRUE, summaryFunction = multiClassSummary,
                                  sampling = subsampleparam  ## adding downsampling here
        )
    }
    
    
    cohort_multi_class_modeltraining <- function(training, testing, outcomelabel, featurelabels, trcontrol, modellabel, seedparam) {
        set.seed(seedparam)
        print(paste0("Training ", modellabel, " model"))
        
        if (modellabel == "xgbTree") {  ## Special params for xgboost
            modfit <- train(x=training[,featurelabels, drop=FALSE], y = training[,outcomelabel], method = modellabel, 
                            # trControl=trcontrol, preProc = c("center", "scale"), metric = "ROC", verbose = FALSE, verbosity = 0)
                            trControl=trcontrol, preProc = c("center", "scale"), metric = "logLoss", verbose = FALSE, verbosity = 0)
        } else if (modellabel == "glmnet") { ## Special trick for lasso with a single variable
            if (ncol(training) == 1) {
                training[,"ones"] <- rep(1, nrow(training))
                testing[,"ones"] <- rep(1, nrow(testing))
                featurelabels <- c(featurelabels, "ones")
            }
            
            
            # A bunch of custom lambda coding:
            ## https://stackoverflow.com/questions/48280074/r-how-to-let-glmnet-select-lambda-while-providing-an-alpha-range-in-caret
            lambda <- unique(unlist(lapply(seq(0, 1, 0.05), function(x){
                init <- glmnet::glmnet(Matrix::as.matrix(training[,c(featurelabels)]), training[,outcomelabel],
                                       family = "multinomial",
                                       # family = "binomial",
                                       nlambda = 100,
                                       alpha = x)
                lambda <- c(min(init$lambda), max(init$lambda))
            }
            )))
            tuneGridb <- expand.grid(.alpha = seq(0, 1, 0.05),
                                     .lambda = seq(min(lambda), max(lambda), length.out = 5))
            
            
            modfit <- suppressWarnings(train(x=training[,c(featurelabels), drop=FALSE], y = training[,outcomelabel], method = modellabel, 
                                             trControl=trcontrol, preProc = c("center", "scale"), metric = "logLoss", verbose = FALSE, verbosity = 0
                                             # tuneGrid = tuneGridb
            ))
        } else if (modellabel == "svmRadial") { ## Theres a bug in svm where you have to use formula instead of x,y ## TOSH WORKING HERE
            modfit <- train(as.formula(paste0(outcomelabel, "~.")), data = training, method = modellabel, 
                            trControl=trcontrol, preProc = c("center", "scale"), metric = "logLoss", trace = FALSE)
        } else {
            modfit <- train(x=training[,featurelabels, drop=FALSE], y = training[,outcomelabel], method = modellabel, 
                            trControl=trcontrol, preProc = c("center", "scale"), metric = "logLoss", trace = FALSE)
        }
        
        ## Train results as well?
        return_model_results <- function(modfit, outcomes_intable, features_intable) {
            predict_labels <- data.frame(predict(modfit, features_intable))
            colnames(predict_labels) <- paste0(modellabel, "_labels_true")
            predict_values <- data.frame(predict(modfit, features_intable, type = "prob"))
            colnames(predict_values) <- c(paste0(outcomelevels, "_pred_", modfit$method))
            # predicted_values_and_labels <- cbind(predict_values, predict_labels)
            
            confstats_roughout <- confusionMatrix(predict_labels[,1], outcomes_intable[,1])
            xx1 <- as.data.frame.matrix(t(confstats_roughout[["table"]]))
            colnames(xx1) <- paste0("predicted_class_", colnames(xx1))
            xx2 <- data.frame(confstats_roughout[["byClass"]])
            rownames(xx2) <- gsub("Class: ", "", rownames(xx2))
            
            # If we only have one class - then we need to dummy and reshape to make it behave
            if (ncol(xx2) == 1){
                xx2 <- rbind(t(xx2), t(xx2))
                rownames(xx2) <- rownames(xx1)
            }
            
            xx3 <- merge(xx1, xx2, by = "row.names")
            rownames(xx3) <- xx3[,"Row.names"]
            xx3 <- xx3[rownames(xx1),]
            confstatsout <- t(xx3[,2:ncol(xx3)])
            confstatsout <- cbind(metric = rownames(confstatsout), confstatsout)
            
            overall_statout <- data.frame(confstats_roughout[[3]])
            colnames(overall_statout) <- modellabel
            overall_statout <- cbind(metric = rownames(overall_statout), overall_statout)
            
            return(list(predict_labels = predict_labels, predict_values = predict_values,
                        confstatsout = confstatsout, overall_statout = overall_statout))
        }
        
        train_model_results <- return_model_results(modfit, outcomes_intable = training[,outcomelabel, drop=FALSE],
                                                    features_intable = training[,featurelabels, drop=FALSE])
        # If we have test results - return, otherwise, return the train results again
        if (nrow(testing[,outcomelabel, drop=FALSE]) > 0){
            test_model_results <- return_model_results(modfit, outcomes_intable = testing[,outcomelabel, drop=FALSE],
                                                       features_intable = testing[,featurelabels, drop=FALSE])
        } else {
            test_model_results <- train_model_results
        }
        
        return(list(train_predict_values_and_labels = cbind(train_model_results$predict_values, train_model_results$predict_labels),
                    train_confstatsout = train_model_results$confstatsout, train_overall_statout = train_model_results$overall_statout, 
                    test_predict_values_and_labels = cbind(test_model_results$predict_values, test_model_results$predict_labels),
                    test_confstatsout = test_model_results$confstatsout, test_overall_statout = test_model_results$overall_statout, 
                    modfit = modfit))
    }
    # train some models
    
    # Im gonna do them separately instead of a for loop for now - probably not as neat, but easier to track for now
    model_out_list <- list()
    if ("glm" %in% models_to_run) {
        glm_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                    modellabel = "glm", seedparam = seedparam) ## generalized linear model
        model_out_list <- c(model_out_list, glm_out = list(glm_out))
    } else { glm_out <- NULL }
    if ("xgbTree" %in% models_to_run) {
        xgbTree_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                        modellabel = "xgbTree", seedparam = seedparam) ## xgboost
        model_out_list <- c(model_out_list, xgbTree_out = list(xgbTree_out))
    } else { xgbTree_out <- NULL }
    if ("svmRadial" %in% models_to_run) {
        svmRadial_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                          modellabel = "svmRadial", seedparam = seedparam) ## SVM with a radial kernel
        model_out_list <- c(model_out_list, svmRadial_out = list(svmRadial_out))
    } else { svmRadial_out <- NULL }
    if ("rf" %in% models_to_run) {
        rf_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                   modellabel = "rf", seedparam = seedparam) ## random forest
        model_out_list <- c(model_out_list, rf_out = list(rf_out))
    } else { rf_out <- NULL }
    if ("glmnet" %in% models_to_run) {
        glmnet_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                       modellabel = "glmnet", seedparam = seedparam) ## lasso
        model_out_list <- c(model_out_list, glmnet_out = list(glmnet_out))
    } else { glmnet_out <- NULL }
    if ("nnet" %in% models_to_run) {
        nnet_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                     modellabel = "nnet", seedparam = seedparam) ## neural net
        model_out_list <- c(model_out_list, nnet_out = list(nnet_out))
    } else { nnet_out <- NULL }
    if ("kknn" %in% models_to_run) {
        kknn_out <- cohort_multi_class_modeltraining(training, testing, outcomelabel, featurelabels, trcontrol,
                                                     modellabel = "kknn", seedparam = seedparam) ## k-nearest neighbor
        model_out_list <- c(model_out_list, kknn_out = list(kknn_out))
    } else { kknn_out <- NULL }
    
    # Failsafe to remove NULL models
    model_out_list <- Filter(Negate(is.null), model_out_list)
    modelnames <- gsub("_out", "", names(model_out_list))
    
    # Collapse train modelstats
    train_modelstats_outtable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "metric", all.x = TRUE, sort = FALSE),
                                        Filter(Negate(is.null), lapply(model_out_list, function(x) x$train_confstatsout)))
    colnames(train_modelstats_outtable) <- c("metric", do.call(paste0, expand.grid(paste0(outcomelevels, "__"), modelnames)))
    train_overallstat_outtable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "metric", all.x = TRUE, sort = FALSE),
                                         Filter(Negate(is.null), lapply(model_out_list, function(x) x$train_overall_statout)))
    
    # Collapse test modelstats
    test_modelstats_outtable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "metric", all.x = TRUE, sort = FALSE),
                                       Filter(Negate(is.null), lapply(model_out_list, function(x) x$test_confstatsout)))
    colnames(test_modelstats_outtable) <- c("metric", do.call(paste0, expand.grid(paste0(outcomelevels, "__"), modelnames)))
    test_overallstat_outtable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "metric", all.x = TRUE, sort = FALSE),
                                        Filter(Negate(is.null), lapply(model_out_list, function(x) x$test_overall_statout)))
    
    # If we have test results - return, otherwise, return the train results again
    if (nrow(testing[,outcomelabel, drop=FALSE]) == 0) {testing <- training}
    
    ## ROC and AUC stats and plots
    ## TRAIN
    ROC_model_data_plot_list <- lapply(model_out_list, function(x) {
        list(ROC_model = x[["modfit"]], ROC_data = training, plot_number = 1, ROC_name = paste0(x[["modfit"]][["method"]], "_train"))
    })
    ROC_custom_out <- ROC_custom_plotter(ROC_model_data_plot_list = ROC_model_data_plot_list, outcome_label = colnames(outcometable), plot_ci = TRUE)
    train_rocplot <- ROC_custom_out[["ROC_plots"]][[1]]
    train_roc_auc_with_ci_res <- cbind(Var = rownames(ROC_custom_out[["ROC_AUC_tables"]][[1]]), ROC_custom_out[["ROC_AUC_tables"]][[1]])
    
    ## TEST
    ROC_model_data_plot_list <- lapply(model_out_list, function(x) {
        list(ROC_model = x[["modfit"]], ROC_data = testing, plot_number = 1, ROC_name = paste0(x[["modfit"]][["method"]], "_test"))
    })
    ROC_custom_out <- ROC_custom_plotter(ROC_model_data_plot_list = ROC_model_data_plot_list, outcome_label = colnames(outcometable), plot_ci = FALSE)
    test_rocplot <- ROC_custom_out[["ROC_plots"]][[1]]
    test_roc_auc_with_ci_res <- cbind(Var = rownames(ROC_custom_out[["ROC_AUC_tables"]][[1]]), ROC_custom_out[["ROC_AUC_tables"]][[1]])
    
    
    # Also want a model list - cause we want to grab the best model to use later on
    modellist <- lapply(model_out_list, function(x) x[["modfit"]])
    names(modellist) <- paste0(modelnames, "_modfit")
    
    # Return variable importance for our models
    varImp_outtable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "metric", all = TRUE, sort = FALSE), 
                              Filter(Negate(is.null), lapply(modellist, function(x) {cbind(metric = rownames(varImp(x)[[1]]), varImp = varImp(x)[[1]][,1,drop=FALSE])})))
    colnames(varImp_outtable) <- c("metric", names(modellist))
    
    return(list(training = training, testing = testing,
                train_modelstats_outtable = train_modelstats_outtable, train_overallstat_outtable = train_overallstat_outtable,
                test_modelstats_outtable = test_modelstats_outtable, test_overallstat_outtable = test_overallstat_outtable,
                train_roc_auc_with_ci_res = train_roc_auc_with_ci_res, train_rocplot=train_rocplot,
                test_roc_auc_with_ci_res = test_roc_auc_with_ci_res, test_rocplot=test_rocplot,
                modellist = modellist, varImp_outtable = varImp_outtable))
    
}



# --------------------------------- helper function to return model results ---------------------------------

return_model_results <- function(modfit, outcomes_intable, features_intable) {
    predict_labels <- data.frame(predict(modfit, features_intable))
    colnames(predict_labels) <- paste0(modellabel, "_labels_true")
    predict_values <- data.frame(predict(modfit, features_intable, type = "prob"))
    colnames(predict_values) <- c(paste0(outcomelevels, "_pred_", modfit$method))
    # predicted_values_and_labels <- cbind(predict_values, predict_labels)
    
    confstats_roughout <- confusionMatrix(predict_labels[,1], outcomes_intable[,1])
    xx1 <- as.data.frame.matrix(t(confstats_roughout[["table"]]))
    colnames(xx1) <- paste0("predicted_class_", colnames(xx1))
    xx2 <- data.frame(confstats_roughout[["byClass"]])
    rownames(xx2) <- gsub("Class: ", "", rownames(xx2))
    
    # If we only have one class - then we need to dummy and reshape to make it behave
    if (ncol(xx2) == 1){
        xx2 <- rbind(t(xx2), t(xx2))
        rownames(xx2) <- rownames(xx1)
    }
    
    xx3 <- merge(xx1, xx2, by = "row.names")
    rownames(xx3) <- xx3[,"Row.names"]
    confstatsout <- t(xx3[,2:ncol(xx3)])
    confstatsout <- cbind(metric = rownames(confstatsout), confstatsout)
    
    overall_statout <- data.frame(confstats_roughout[[3]])
    colnames(overall_statout) <- modellabel
    overall_statout <- cbind(metric = rownames(overall_statout), overall_statout)
    
    return(list(predict_labels = predict_labels, predict_values = predict_values,
                confstatsout = confstatsout, overall_statout = overall_statout))
}









## Return prediction stats from confmat with a vector of true labels and test labels
return_pred_stats <- function(true_labels, test_labels) {
    confstats_roughout <- confusionMatrix(true_labels, test_labels)
    confstats_tempttab <- merge(as.data.frame.matrix(t(confstats_roughout[["table"]])), 
                                data.frame(confstats_roughout[["byClass"]], row.names = gsub("Class: ", "", rownames(data.frame(confstats_roughout[["byClass"]])))), 
                                by = "row.names")
    confstatsout <- t(data.frame(confstats_tempttab[,2:ncol(confstats_tempttab)], row.names = paste0("predicted_class_", colnames(confstats_roughout[["table"]]))))
    confstatsout <- cbind(metric = rownames(confstatsout), confstatsout)
    overall_statout <- data.frame(confstats_roughout[[3]])
    overall_statout <- cbind(metric = rownames(data.frame(confstats_roughout[[3]])), data.frame(confmat_outstat = confstats_roughout[[3]]))
    return(list(confstatsout = confstatsout, overall_statout = overall_statout))
}





# --------------------------------- CUSTOM ROC PLOTTER ---------------------------------
# ROC_model_data_plot_list <- list(
#     combo1 = list(ROC_model = mn.net, ROC_data = iris_threeclass_test, plot_number = 1, ROC_name = "mn_test"),
#     combo2 = list(ROC_model = mn.net, ROC_data = iris_threeclass_train, plot_number = 1, ROC_name = "mn_train"),
#     combo3 = list(ROC_model = rfmodel, ROC_data = iris_threeclass_test, plot_number = 2, ROC_name = "rf_test"),
#     combo4 = list(ROC_model = rfmodel, ROC_data = iris_threeclass_train, plot_number = 2, ROC_name = "rf_train")
# )
# outcome_label = "Species"
ROC_custom_plotter <- function(ROC_model_data_plot_list, outcome_label, plot_ci = TRUE) {
    # multiclass = FALSE ## This can be detected automatically
    
    # First grab all of the ROCs we have from our model data combos
    ROC_outlist_all <- list()
    multiclass_AUC_list <- list()
    for (model_data_combo in seq_len(length(ROC_model_data_plot_list))) {
        grab_model <- ROC_model_data_plot_list[[model_data_combo]][["ROC_model"]]
        grab_data <- ROC_model_data_plot_list[[model_data_combo]][["ROC_data"]]
        ROC_predictions <- predict(grab_model, grab_data, type = "prob")
        ROC_values_out <- multiclass.roc(grab_data[,outcome_label], ROC_predictions)
        ROC_objects <- ROC_values_out[["rocs"]]
        if (length(ROC_objects) == 1) { # Then its a two class
            ## Ok - so if its two class by definition - then even though some models will output both sides A > B and B > A
            ## We only need the first regardless, so lets just double check that we grab the first one always?
            grab_ROC <- ROC_objects[[1]]
            if (length(grab_ROC) == 2) {grab_ROC <- grab_ROC[[1]]}
            ROC_outlist_all[[paste0("combo", model_data_combo)]] <- list(ROC = grab_ROC,
                                                                         plot_number = ROC_model_data_plot_list[[model_data_combo]][["plot_number"]],
                                                                         ROC_name = ROC_model_data_plot_list[[model_data_combo]][["ROC_name"]])
        } else { # Multiclass!
            grab_ROC <- lapply(ROC_objects, function(x) {
                out1 <- x[1]
                names(out1) <- "ROC"
                subplotname <- paste(x[[1]][["levels"]], collapse = "/")
                out2 <- list(plot_number = ROC_model_data_plot_list[[model_data_combo]][["plot_number"]],
                             ROC_name = paste0(ROC_model_data_plot_list[[model_data_combo]][["ROC_name"]], "_", subplotname))
                c(out1, out2)
            }) ## Just grab one of the two sides for each one # TOSH WORKING HERE
            ## so i guess we still put these all in one plot - and up to the user to not put them all in one place - cause that would be dumb (or maybe not)
            ROC_outlist_all <- c(ROC_outlist_all, grab_ROC)
            multiclass_AUC <- list(ROC_name = ROC_model_data_plot_list[[model_data_combo]][["ROC_name"]],
                                   plot_number = ROC_model_data_plot_list[[model_data_combo]][["plot_number"]],
                                   multiclass_AUC = ROC_values_out[["auc"]][[1]])
            multiclass_AUC_list[[model_data_combo]] <- multiclass_AUC
            ## Ok - but we need to store the multiclass AUC if applicable
            
        }
    }
    multiclass_AUC_table <- do.call(rbind, multiclass_AUC_list)
    ## Now that we have all of our ROCs, lets split this into groups by plotnumber
    plotgroups <- paste0("plot", unlist(lapply(ROC_outlist_all, function(x) x$plot_number)))
    # ROC_outlist_byplot <- lapply(split(ROC_outlist_all, f = plotgroups), function(x) 
    ROC_outlist_byplot <- lapply(split(ROC_outlist_all, f = plotgroups), function(x) {
        out1 <- lapply(x, function(y) {y[["ROC"]]})
        names(out1) <- lapply(x, function(y) {y[["ROC_name"]]})
        out1
    })
    # ROC_outlist_byplot
    # ROC_outlist_byplot[["plot1"]]
    
    # Now for each plot - do the analyses we want to do:
    outplot_list <- list()
    ROC_AUC_list <- list()
    for (plotnumber_list_num in seq_len(length(ROC_outlist_byplot))) {
        
        # Grab the ROC list for this plot
        grab_ROC_list_object <- ROC_outlist_byplot[[plotnumber_list_num]]
        
        # Get AUC info
        ROC_AUC_vals <- data.frame(lapply(grab_ROC_list_object, function(x) auc(x)[[1]]), row.names = "AUC_val")
        ROC_AUC_CI_vals <- data.frame(do.call(cbind, suppressWarnings(lapply(grab_ROC_list_object, function(x) ci.auc(x)[c(1,3)]))), 
                                      row.names = c("lowerCI", "upperCI"))
        ROC_AUC_outtable <- data.frame(t(rbind(ROC_AUC_vals, ROC_AUC_CI_vals)))
        
        # Check for a multiclass AUC - and if so then insert here
        if (!is.null(multiclass_AUC_table) & sum(multiclass_AUC_table[,"plot_number"] %in% plotnumber_list_num) > 0){
            multiclass_AUC_table_sel <- multiclass_AUC_table[multiclass_AUC_table[,"plot_number"] %in% plotnumber_list_num,,drop=FALSE]
            ROC_AUC_outtable[,"multiclass_AUC"] <- unlist(apply(multiclass_AUC_table_sel, 1, function(x)
                rep(x[3], sum(grepl(x[1], rownames(ROC_AUC_outtable))))))
        }
        ROC_AUC_list[[plotnumber_list_num]] <- ROC_AUC_outtable
        
        # Create AUC labels to attach to the names for plotting
        model_AUC_labels <- apply(ROC_AUC_outtable, 1, function(x) {
            paste0("\nAUC: ", round(x[1], 3), " (", round(x[2], 3), "-", round(x[3], 3), ")")
        })
        names(grab_ROC_list_object) <- paste0(names(model_AUC_labels), model_AUC_labels)
        
        # Make our simple plot
        pout <- ggroc(grab_ROC_list_object)
        pout <- pout + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey")
        pout <- pout + theme_minimal()  + coord_equal()
        
        # Plot CI if param is set to yes
        if (plot_ci) {
            # First calculate the cis for it
            ci.list <- lapply(grab_ROC_list_object, function(x) ci.se(x, specificities = seq(0, 1, l = 25), boot.n = 1000))
            dat.ci.list <- lapply(ci.list, function(ciobj) {data.frame(x = as.numeric(rownames(ciobj)), lower = ciobj[, 1], upper = ciobj[, 3])})
            # Add CIs
            for (model_ci_num in seq_len(length(dat.ci.list))) {
                pout <- pout + geom_ribbon(
                    data = dat.ci.list[[model_ci_num]],
                    aes(x = x, ymin = lower, ymax = upper),
                    fill = model_ci_num + 1, alpha = 0.2, inherit.aes = F) 
            }
        }
        outplot_list[[plotnumber_list_num]] <- pout
    }
    
    return(list(ROC_plots = outplot_list, ROC_AUC_tables = ROC_AUC_list))
}



# --------------------------------- Combine the AUC and stattable ---------------------------------
combine_mgc_modelstats_and_auc_table <- function(modelstats_outtable, roc_auc_with_ci_res) {
    modelstats_clean <- data.frame(modelstats_outtable[,2:ncol(modelstats_outtable)], row.names = modelstats_outtable[,"metric"])
    colnames(modelstats_clean) <- paste0(c("c1_", "c2_"), unlist(lapply(strsplit(colnames(modelstats_clean), split = "__"), 
                                                                              function(x) x[length(x)])))
    auc_stats_clean_temp <- data.frame(t(data.frame(roc_auc_with_ci_res, 
        row.names = unlist(lapply(strsplit(rownames(roc_auc_with_ci_res), split = "_"), function(x) x[1])))[,c(2:4)]), row.names = c("AUC", "AUClowerCI", "AUCupperCI"))
    auc_stats_clean <- auc_stats_clean_temp[,rep(1:ncol(auc_stats_clean_temp), each = 2)]
    colnames(auc_stats_clean) <- apply(expand.grid(c("c1_", "c2_"), colnames(auc_stats_clean_temp)), 1, paste, collapse="")
    allstat_table <- rbind(modelstats_clean, auc_stats_clean)
    allstat_table[] <- apply(allstat_table, 2, as.numeric)
    return(allstat_table)
}



# --------------------------------- END ---------------------------------





