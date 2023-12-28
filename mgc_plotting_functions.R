################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Plotting Functions

## Load in Libraries
packagelist = c("ggplot2", "reshape2", "DESeq2", "grid", "gridExtra", "scales", "ggrepel", "circlize", "viridis", "RColorBrewer", "tools", "ggpubr", "ggsignif", "ComplexHeatmap", "stringr", "tidyverse", "ggsci", "Rtsne", "GGally", "corrplot", "waffle", "ggthemes")
junk <- lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
ht_opt$message = FALSE

# https://nanx.me/ggsci/reference/pal_npg.html
## npg_color_list: c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")

# --------------------------------- barchart ---------------------------------

# infile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/newman-pace/output_newcohorts_age_covar/gsea/comp_cohort__pace_v_healthy/comp_cohort__pace_v_healthy_gsea_HALL.csv"
# intable <- read.table(infile, sep = ",", header = TRUE)
# bartab <- intable[intable[,"p.adjust"] < 0.1,c(1,6,8)]
# labsparam <- list(title = "top pathways", x = "pathway", y = "NES", fill = "Adjpvalue")
# pout <- plot_barchart(bartab, colorcol = 3, labsparam = labsparam)
# errorbarparam - what to do error bars off of? {ci, sd (standard dev)} ## ADD MORE OPTIONS AS NEEDED

plot_barchart <- function(bartab, colorvar = NULL, errorbarvar = NULL, labsparam, strwrap_characternum = 20){

    plotdata <- bartab
    ## First - replace the label strings with spaces so we can wrap nicely
    plotdata[,1] <- gsub("_", " ", plotdata[,1])
    ## Stop reordering
    plotdata[,1] <- factor(plotdata[,1], levels = plotdata[,1])

    ## Add colors if input
    if (!is.null(colorvar)) {plotdata$colorvar = colorvar[,1]}

    pout <- ggplot(plotdata, mapping = aes(x = str_wrap(plotdata[,1], strwrap_characternum), y = plotdata[,2]))
    # if (!is.null(colorcol)) {pout <- ggplot(plotdata, mapping = aes(x = str_wrap(plotdata[,1], strwrap_characternum), y = plotdata[,2], fill = plotdata[,colorcol]))}
    if (!is.null(colorvar)) {
        if (is.numeric(plotdata[,3])){
            pout <- pout + aes(fill=plotdata$colorvar)
            colormin <- min(plotdata$colorvar, 0)
            colormax <- max(plotdata$colorvar, 0)
            pout <- pout + scale_fill_viridis(limits = c(colormin, colormax), direction = -1, option = "viridis")
        } else {
            pout <- pout + aes(fill=plotdata$colorvar)
            ## If literally a color - then use the color - 2022-11-03
            if (sum(areColors(plotdata$colorvar)) == length(plotdata$colorvar) & !is.numeric(plotdata$colorvar)) { ## If your colorvar are literally colors and NOT number - then use them
                pout <- pout + scale_color_identity()
                pout <- pout + scale_fill_identity()
            }
            ## If literally a color - then use the color - 2022-11-03
        }
    }
    pout <- pout + geom_bar(stat = "identity")
    # pout <- pout + scale_x_discrete(limits=rev(str_wrap(plotdata[,1], strwrap_characternum)))
    pout <- pout + scale_x_discrete(limits=str_wrap(plotdata[,1], strwrap_characternum))
    if (!is.null(errorbarvar)) {
        pout <- pout + geom_errorbar(aes(ymax = plotdata[,2] + errorbarvar[,1], ymin = plotdata[,2] - errorbarvar[,1]),
                                     width = 0.5, position = position_dodge(width = 0.9))}
    pout <- pout + labs(x = labsparam$x, y = labsparam$y, fill = labsparam$fill, title = labsparam$title)
    pout <- pout + theme_pubr(x.text.angle = 90, legend = "right")
    pout
    return(pout)
}



# Make the table for a boxplot
# indata - standard data table, samples on y, features on x, with row and col names
# inmeta - standard meta data, samples on x, categories on y, col names, FIRST COLUMN IS SAMPLES
# FOI - features of interest, with values to be plotted in the boxplot
# COI - category of interest - the primary category to split data on
# COIaddon - categorical information that will be added onto features if you want subsections of features
#   ex) C1F1, C1F2, C2F1, C2F2 for 4 boxplots
make_boxplot_table <- function(indata, inmeta, FOI, COI, COIaddon=NULL) {
    ## Combine the data and metadata into one table
    temp = rbind(t(inmeta[match(inmeta[,1], colnames(indata)),c(2:ncol(inmeta)), drop=FALSE]), as.matrix(indata))
    colnames(temp) = colnames(indata)
    # tempmelt = as.matrix(melt(temp, varnames = c("vars","sample")))
    tempmelt = melt(temp, varnames = c("vars","sample"))

    ## Select the primary category of interest and replace samples with it
    tempmelt[,2] = inmeta[,COI][match(tempmelt[,2], inmeta[,1])]
    ## Add on any data from COIaddon to features
    if (!is.null(COIaddon)) {tempmelt[,1] =
        paste(as.character(tempmelt[tempmelt[,1]==COIaddon,3]), as.character.factor(tempmelt[,1]), sep = "_")}

    ## Select the FOI
    featureindex = c()
    for (FOIcount in seq_len(length(FOI))){
        # featureindex = c(featureindex, grep(FOI[FOIcount], tempmelt[,1], fixed = TRUE))
        ## Return the index that has the GOI from the melted table
        # featureindex = c(featureindex, rownames(tempmelt[tempmelt[,1] == FOI,])) # edit 20200410
        featureindex = c(featureindex, rownames(tempmelt[tempmelt[,1] == FOI[FOIcount],]))

    }
    tempmelt = tempmelt[featureindex,]
    tempmelt[,1] = as.character(tempmelt[,1])
    tempmelt[,2] = as.character(tempmelt[,2])
    tempmelt[,3] = as.numeric(as.character(tempmelt[,3]))
    return(data.frame(tempmelt))
}


# --------------------------------- monster function for boxplots ---------------------------------

## GGplot the boxplots for SOI for each gene in the item - MAJOR REWRITE - but necessary.
# Input the alread melted boxplot (function above)
# if want to label with second axis, then NAME the second axis, and the second axis will be scaled to the first
#   ex) group1 and group2, if you NAME group2, then group1 will get the MAJOR markings, and group2 gets the SECOND markings
# xsplit - is the param by which you split your x-axis on <"category" | "feature" > Do you want a boxplot(s) for each category or for each feature? feature is the feature in your count table, category is the grouping in your metadata
# colorparam - named vector with the names being the features, and the named item being the color you want it
# plotstats - < "inter", "intra" > to say if you want stats for INTERgroup comparisons or INTRAgroup comparisons - intra will automatically facet the plot, inter will leave as one plot
# comparisonparam - if inputted, should be a list of features/categories to compare
###### TOSH - IT GOES FEATURE, CATEGORY, VALUE. this may not be how you want it, but it is what it is.
# labsparam - title, x, y, fill, catorder, featorder
# boxplot_plotter <- function(boxplottable,
#                             xsplit = "category",
#                             labsparam, colorparam = NULL, secondaxis=NULL,
#                             plotstats = FALSE, comparisonparam = NULL, testtypeparam = "wilcox.test", pointcoloring = FALSE, 
#                             plotpoints = TRUE ## Added 2022-08-01
#                             ) {
# 
#     ## INPUTTING THE PREMADE MELTED TABLE
#     # ## Turn counttab and and metatable into melted plottab
#     # ## Melt the table down
#     # meltedtab = data.frame(melt(as.matrix(intable), varnames = c("features", "samples")), stringsAsFactors = FALSE)
#     # ## Reassign the samples as the class
#     # meltedtab[,2] = metatable[,2][match(meltedtab[,2], metatable[,1])]
#     # meltedtab[,1] = as.character(meltedtab[,1])
#     # plottab = meltedtab
#     plottab = boxplottable
# 
#     ## Create failsafe if you omit the catorder and featorder ##28-07-2021
#     if (is.null(labsparam$featorder)) {labsparam$featorder <- unique(plottab[,1])}
#     if (is.null(labsparam$catorder)) {labsparam$catorder <- unique(plottab[,2])}
#     
#     ## To get numbers in the same world for plotting purpose - we are going to round them to the nearest base 10
#     ## Really only usable when you have two categories - PROBABLY BROKEN UNTIL FURTHER TESTED
#     ## use the NAMED secondaxis variable, and scale that variable up to the first one, while labeling the second axis
#     if (!is.null(secondaxis)) {
#         roundUp <- function(x) 10^ceiling(log10(x))
#         for (i in seq_len(length(secondaxis))) {
#             temptab1 = plottab[plottab[,1] == secondaxis[i],]
#             temptab2 = plottab[plottab[,1] != secondaxis[i],]
#             tempval = mean(na.omit(temptab2[,3])) / mean(na.omit(temptab1[,3]))
#             scalefactor = roundUp(tempval)
#             plottab[plottab[,1] == secondaxis[i],3] = plottab[plottab[,1] == secondaxis[i],3] * scalefactor
#         }
#     }
# 
#     ## Get maxval on plot - to be used for plotting purposes later
#     maxval = round(max(na.omit(plottab[,3])),2)
#     minval = round(min(na.omit(plottab[,3])),2)
#     valrange <- abs(maxval) + abs(minval)
# 
#     if (xsplit == "feature") {
#         xcol = 1
#         zcol = 2
#         ## 3/17/12
#         # bporder = labsparam$featorder
#         # bporderintra = labsparam$catorder
#         bporder = factor(labsparam$featorder, levels = labsparam$featorder)
#         bporderintra = factor(labsparam$catorder, levels = labsparam$catorder)
#         plottab[,xcol] <- factor(plottab[,xcol], levels =  labsparam$featorder)
#         plottab[,zcol] <- factor(plottab[,zcol], levels =  labsparam$catorder)
#         ## 3/17/12
#         if(is.null(bporder)) {bporder <- unique(plottab[,zcol])}
#         if(is.null(bporderintra)) {bporderintra <- unique(plottab[,xcol])}
#     } else { # xsplit == "category"
#         xcol = 2
#         zcol = 1
#         ## 3/17/12
#         # bporder = labsparam$catorder
#         # bporderintra = labsparam$featorder 
#         bporder = factor(labsparam$catorder, levels = labsparam$catorder)
#         bporderintra = factor(labsparam$featorder, levels = labsparam$featorder)
#         plottab[,xcol] <- factor(plottab[,xcol], levels =  labsparam$catorder)
#         plottab[,zcol] <- factor(plottab[,zcol], levels =  labsparam$featorder)
#         ## 3/17/12
#         if(is.null(bporder)) {bporder <- unique(plottab[,xcol])}
#         if(is.null(bporderintra)) {bporderintra <- unique(plottab[,zcol])}
#     }
#     
#     ## Add in point coloring if desired # added 9/21/2021
#     if (pointcoloring == TRUE) {
#       # pointcolorcol <- 4
#       pointcolorparam <- plottab[,4]
#     } else {
#       # pointcolorcol <- zcol
#       pointcolorparam <- "black"
#     }
#     ## STANDARD BOXPLOT WITH NO STATS
#     pout <- ggplot(plottab, aes(x=plottab[,xcol], y=plottab[,3], z=plottab[,zcol], fill = plottab[,zcol]))
#     pout <- pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
#     # pout <- pout + geom_point(data = plottab, aes(x=plottab[,xcol], y=plottab[,3], fill = plottab[,pointcolorcol]),
#     #                           position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.8, shape = 16)
#     if (plotpoints) {
#         pout <- pout + geom_point(data = plottab, aes(x=plottab[,xcol], y=plottab[,3], color = pointcolorparam),
#                                   position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.8, shape = 16)
#     }
#     
# 
#     pout <- pout + scale_x_discrete(limits = bporder)
#     # pout <- pout + scale_fill_discrete(limits = bporderintra)
# 
#     pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, fill = labsparam$fill)
#     if(!is.null(secondaxis)) {pout <- pout +
#         scale_y_continuous(sec.axis = sec_axis(~./scalefactor, name = secondaxis[1]))}
#     if(!is.null(colorparam)) {pout <- pout + scale_fill_manual(values = unname(colorparam), limits = names(colorparam))}
#     pout <- pout + theme_pubr(x.text.angle = 45, base_size = 12)
# 
#     pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*1.2, digits = 2)),
#                                       # breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*1.2, digits = 2),length.out = 5))
#                                       breaks = round(seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*1.2, digits = 2),length.out = 5),2))
# 
# 
#     ## CREATE COMBO TABLE FOR STATS IF APPLICABLE
#     if (plotstats != FALSE) {
#         if (is.null(comparisonparam)) {
#             combcol = ifelse(plotstats == "inter", xcol, zcol)
#             ## Put a failsafe in here to only test against those with at least 2 points (cant have a t.test against 1 pt)
#             # comparisontemp <- na.omit(names(table(plottab[!is.na(plottab[,combcol]),combcol])[table(plottab[!is.na(plottab[,combcol]),combcol])> 1]))
#             # comparisontemp2 <- comparisontemp[order(comparisontemp)]
# 
#             # if (length(comparisontemp2) > 1) {
#             #     combtab = combn(comparisontemp2, 2, simplify=F)
#             # } else {
#             #     combtab <- NULL
#             # }
#             
#             combtab = combn(as.character(unique(plottab[!is.na(plottab[,combcol]),combcol])), 2, simplify=F)
#             
#         } else {
#             if (is.character(comparisonparam) | is.factor(comparisonparam)) {
#                 comparisonparam <- as.character(comparisonparam)
#                 combtab = combn(comparisonparam, 2, simplify=F)
#             }
#             if (is.list(comparisonparam)) {
#                 combtab = comparisonparam
#             }
#         }
#         ## Failsafe here to exclude statistical testing if not enough values in group
#         excludegroups <- names(table(plottab[,2])[table(plottab[,2]) < 2])
#         combtab <- combtab[!unlist(lapply(combtab, function(x) sum(x %in% excludegroups) > 0))]
#     }
# 
#     ### THIS DOES SIGNIFICANCE BETWEEN THE X AXIS GROUPS
#     if (plotstats == "inter"){
#         # pout <- pout + geom_signif(comparisons = combtab, step_increase = 0.1,
#         #                     y_position = seq(maxval*1.2, (maxval*(1.1+(0.1*length(combtab)))), 0.1*maxval),
#         #                     map_signif_level=FALSE)}
#         pout <- pout + geom_signif(comparisons = combtab, step_increase = 0.15, test = testtypeparam,
#                                    # y_position = seq(maxval*1.1, (maxval*(1.0+(0.1*length(combtab)))), 0.1*maxval), ## 27-07-2021 - change because the step.increase does this better than the manual setting
#                                    map_signif_level=FALSE)
#         pout <- suppressMessages(pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.15*length(combtab))), digits = 2)),
#                                       # breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5)))
#                                       breaks = round(seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5),2)))
#     }
# 
# 
#     ### THIS DOES SIGNIFICANCE BETWEEN SUBGROUPS
#     if (plotstats == "intra") {
#         pout <- ggplot(plottab, aes(x=plottab[,zcol], y=plottab[,3], fill = plottab[,zcol]))
#         pout <- pout + facet_wrap(~factor(plottab[,xcol], levels = bporder), scales = "free")
#         pout <- pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
#         # pout <- pout + geom_point(data = plottab, aes(x=plottab[,zcol], y=plottab[,3], fill = plottab[,pointcolorcol]),
#         #                           position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.5, shape = 16)
#         if (plotpoints) {
#             pout <- pout + geom_point(data = plottab, aes(x=plottab[,zcol], y=plottab[,3], color = pointcolorparam),
#                                       position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.5, shape = 16)
#         }
#         
# 
# 
#         pout <- pout + scale_x_discrete(limits = bporderintra)
# 
#         pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, fill = labsparam$fill)
#         if(!is.null(secondaxis)) {pout <- pout +
#             scale_y_continuous(sec.axis = sec_axis(~./scalefactor, name = secondaxis[1]))}
#         if(!is.null(colorparam)) {pout <- pout + scale_fill_manual(values = unname(colorparam), limits = names(colorparam))}
#         pout <- pout + theme_pubr(x.text.angle = 60, base_size = 12, border = TRUE)
#         # pout <- pout + stat_compare_means(comparisons = combtab, method.args = list(exact = FALSE), label = "p.format", size = 3,
#         #                        label.y = seq(maxval*1.1, (maxval*(1.0+(0.1*length(combtab)))), 0.1*maxval))
#         if (length(combtab) > 0) { # failsafe for if there arent actually nay viable comps
#             ## Failsafe for single comp and fix rounding error:
#             failsafesignparam <- sign(maxval+(0.1*length(combtab)*valrange) - maxval*1.1)
#             pout <- pout + stat_compare_means(comparisons = combtab, method = testtypeparam,
#                                               # method.args = list(exact = TRUE), ### FOR SOME REASON - THIS IS KILLING MY RUN TIME, so omitting for now
#                                               label = "p.format", size = 3,
#                                               # label.y = seq(maxval*1.1, (maxval*(1.0+(0.1*length(combtab)))), 0.1*maxval))
#                                               # label.y = seq(maxval*1.1, maxval*(0.1*length(combtab)*valrange), 0.1*valrange))
#                                               label.y = seq(maxval*1.1, maxval+(0.1*length(combtab)*valrange), failsafesignparam*0.1*valrange))
#                                               ## bug here where the rounding screws it up... really only when you have one comp though
#                                               # label.y = seq(ceiling(maxval*1.1, 2), ceiling(maxval+(0.1*length(combtab)*valrange),2), round(0.1*valrange,2)))
#         }
# 
#         # pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2)),
#         #                     breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5))
#         # pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(0.1*length(combtab)*valrange, digits = 2)),
#         #                                   breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5))
#         pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval+(0.1*length(combtab)*valrange+0.1*valrange), digits = 2)),
#                                           # breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5))
#                                           breaks = round(seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5),2))
# 
#         # expand_scale(mult = c(1.2, round(maxval*(1.1+(0.1*length(combtab))), digits = 2)))
# 
#     }
# 
#     ## Add in point coloring if desired # added 9/21/2021
#     if (pointcoloring != TRUE) {
#       pout <- pout + scale_color_manual(values = "black", breaks = "black")
#     }
#     
#     pout
# 
#     #return(pout)
# 
# }

boxplot_plotter <- function(boxplottable,
                            xsplit = "category",
                            labsparam, colorparam = NULL, secondaxis=NULL,
                            plotstats = FALSE, comparisonparam = NULL, testtypeparam = "wilcox.test", pointcoloring = FALSE, 
                            plotpoints = TRUE, ## Added 2022-08-01
                            violin = FALSE ## Added 2022-10-17
) {
    
    ## INPUTTING THE PREMADE MELTED TABLE
    # ## Turn counttab and and metatable into melted plottab
    # ## Melt the table down
    # meltedtab = data.frame(melt(as.matrix(intable), varnames = c("features", "samples")), stringsAsFactors = FALSE)
    # ## Reassign the samples as the class
    # meltedtab[,2] = metatable[,2][match(meltedtab[,2], metatable[,1])]
    # meltedtab[,1] = as.character(meltedtab[,1])
    # plottab = meltedtab
    plottab = boxplottable
    
    ## Create failsafe if you omit the catorder and featorder ##28-07-2021
    if (is.null(labsparam$featorder)) {labsparam$featorder <- unique(plottab[,1])}
    if (is.null(labsparam$catorder)) {labsparam$catorder <- unique(plottab[,2])}
    
    ## To get numbers in the same world for plotting purpose - we are going to round them to the nearest base 10
    ## Really only usable when you have two categories - PROBABLY BROKEN UNTIL FURTHER TESTED
    ## use the NAMED secondaxis variable, and scale that variable up to the first one, while labeling the second axis
    if (!is.null(secondaxis)) {
        roundUp <- function(x) 10^ceiling(log10(x))
        for (i in seq_len(length(secondaxis))) {
            temptab1 = plottab[plottab[,1] == secondaxis[i],]
            temptab2 = plottab[plottab[,1] != secondaxis[i],]
            tempval = mean(na.omit(temptab2[,3])) / mean(na.omit(temptab1[,3]))
            scalefactor = roundUp(tempval)
            plottab[plottab[,1] == secondaxis[i],3] = plottab[plottab[,1] == secondaxis[i],3] * scalefactor
        }
    }
    
    ## Get maxval on plot - to be used for plotting purposes later
    maxval = round(max(na.omit(plottab[,3])),2)
    minval = round(min(na.omit(plottab[,3])),2)
    valrange <- abs(maxval) + abs(minval)
    
    if (xsplit == "feature") {
        xcol = 1
        zcol = 2
        ## 3/17/12
        # bporder = labsparam$featorder
        # bporderintra = labsparam$catorder
        bporder = factor(labsparam$featorder, levels = labsparam$featorder)
        bporderintra = factor(labsparam$catorder, levels = labsparam$catorder)
        plottab[,xcol] <- factor(plottab[,xcol], levels =  labsparam$featorder)
        plottab[,zcol] <- factor(plottab[,zcol], levels =  labsparam$catorder)
        ## 3/17/12
        if(is.null(bporder)) {bporder <- unique(plottab[,zcol])}
        if(is.null(bporderintra)) {bporderintra <- unique(plottab[,xcol])}
    } else { # xsplit == "category"
        xcol = 2
        zcol = 1
        ## 3/17/12
        # bporder = labsparam$catorder
        # bporderintra = labsparam$featorder 
        bporder = factor(labsparam$catorder, levels = labsparam$catorder)
        bporderintra = factor(labsparam$featorder, levels = labsparam$featorder)
        plottab[,xcol] <- factor(plottab[,xcol], levels =  labsparam$catorder)
        plottab[,zcol] <- factor(plottab[,zcol], levels =  labsparam$featorder)
        ## 3/17/12
        if(is.null(bporder)) {bporder <- unique(plottab[,xcol])}
        if(is.null(bporderintra)) {bporderintra <- unique(plottab[,zcol])}
    }
    
    ## Add in point coloring if desired # added 9/21/2021
    if (pointcoloring == TRUE) {
        # pointcolorcol <- 4
        pointcolorparam <- plottab[,4]
    } else {
        # pointcolorcol <- zcol
        pointcolorparam <- "black"
    }
    ## STANDARD BOXPLOT WITH NO STATS
    pout <- ggplot(plottab, aes(x=plottab[,xcol], y=plottab[,3], z=plottab[,zcol], fill = plottab[,zcol]))
    if (!violin) {
        pout <- pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
    } else {
        pout <- pout + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE)
    }
    
    # pout <- pout + geom_point(data = plottab, aes(x=plottab[,xcol], y=plottab[,3], fill = plottab[,pointcolorcol]),
    #                           position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.8, shape = 16)
    if (plotpoints) {
        pout <- pout + geom_point(data = plottab, aes(x=plottab[,xcol], y=plottab[,3], color = pointcolorparam),
                                  position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.8, shape = 16)
    }
    
    
    pout <- pout + scale_x_discrete(limits = bporder)
    # pout <- pout + scale_fill_discrete(limits = bporderintra)
    
    pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, fill = labsparam$fill)
    if(!is.null(secondaxis)) {pout <- pout +
        scale_y_continuous(sec.axis = sec_axis(~./scalefactor, name = secondaxis[1]))}
    if(!is.null(colorparam)) {pout <- pout + scale_fill_manual(values = unname(colorparam), limits = names(colorparam))}
    pout <- pout + theme_pubr(x.text.angle = 45, base_size = 12)
    
    ## 2022-12-13 COMBO OF THESE THINGS is fucking this up again, gah, so... maybe we dont need to scale of the geom_signif works well????
    if (plotstats == FALSE) {
        pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*1.2, digits = 2)),
                                          # breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*1.2, digits = 2),length.out = 5))
                                          breaks = round(seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*1.2, digits = 2),length.out = 5),2))
    }
    
    ## CREATE COMBO TABLE FOR STATS IF APPLICABLE
    if (plotstats != FALSE) {
        if (is.null(comparisonparam)) {
            combcol = ifelse(plotstats == "inter", xcol, zcol)
            ## Put a failsafe in here to only test against those with at least 2 points (cant have a t.test against 1 pt)
            # comparisontemp <- na.omit(names(table(plottab[!is.na(plottab[,combcol]),combcol])[table(plottab[!is.na(plottab[,combcol]),combcol])> 1]))
            # comparisontemp2 <- comparisontemp[order(comparisontemp)]
            
            # if (length(comparisontemp2) > 1) {
            #     combtab = combn(comparisontemp2, 2, simplify=F)
            # } else {
            #     combtab <- NULL
            # }
            
            combtab = combn(as.character(unique(plottab[!is.na(plottab[,combcol]),combcol])), 2, simplify=F)
            
        } else {
            if (is.character(comparisonparam) | is.factor(comparisonparam)) {
                comparisonparam <- as.character(comparisonparam)
                combtab = combn(comparisonparam, 2, simplify=F)
            }
            if (is.list(comparisonparam)) {
                combtab = comparisonparam
            }
        }
        ## Failsafe here to exclude statistical testing if not enough values in group
        excludegroups <- names(table(plottab[,2])[table(plottab[,2]) < 2])
        combtab <- combtab[!unlist(lapply(combtab, function(x) sum(x %in% excludegroups) > 0))]
    }
    
    ### THIS DOES SIGNIFICANCE BETWEEN THE X AXIS GROUPS
    if (plotstats == "inter"){
        # pout <- pout + geom_signif(comparisons = combtab, step_increase = 0.1,
        #                     y_position = seq(maxval*1.2, (maxval*(1.1+(0.1*length(combtab)))), 0.1*maxval),
        #                     map_signif_level=FALSE)}
        
        
        pout <- pout + geom_signif(comparisons = combtab, step_increase = 0.15, test = testtypeparam,
                                   # y_position = seq(maxval*1.1, (maxval*(1.0+(0.1*length(combtab)))), 0.1*maxval), ## 27-07-2021 - change because the step.increase does this better than the manual setting
                                   map_signif_level=FALSE)
        ## 2022-12-13 COMBO OF THESE THINGS is fucking this up again, gah, so... maybe we dont need to scale of the geom_signif works well????
        # pout <- suppressMessages(pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.15*length(combtab))), digits = 2)),
        #                                                    # breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5)))
        #                                                    breaks = round(seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5),2)))
    }
    
    
    ### THIS DOES SIGNIFICANCE BETWEEN SUBGROUPS
    if (plotstats == "intra") {
        pout <- ggplot(plottab, aes(x=plottab[,zcol], y=plottab[,3], fill = plottab[,zcol]))
        pout <- pout + facet_wrap(~factor(plottab[,xcol], levels = bporder), scales = "free")
        if (!violin) {
            pout <- pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
        } else {
            pout <- pout + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE)
        }
        # pout <- pout + geom_point(data = plottab, aes(x=plottab[,zcol], y=plottab[,3], fill = plottab[,pointcolorcol]),
        #                           position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.5, shape = 16)
        if (plotpoints) {
            pout <- pout + geom_point(data = plottab, aes(x=plottab[,zcol], y=plottab[,3], color = pointcolorparam),
                                      position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 0.5, shape = 16)
        }
        
        
        
        pout <- pout + scale_x_discrete(limits = bporderintra)
        
        pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, fill = labsparam$fill)
        if(!is.null(secondaxis)) {pout <- pout +
            scale_y_continuous(sec.axis = sec_axis(~./scalefactor, name = secondaxis[1]))}
        if(!is.null(colorparam)) {pout <- pout + scale_fill_manual(values = unname(colorparam), limits = names(colorparam))}
        pout <- pout + theme_pubr(x.text.angle = 60, base_size = 12, border = TRUE)
        # pout <- pout + stat_compare_means(comparisons = combtab, method.args = list(exact = FALSE), label = "p.format", size = 3,
        #                        label.y = seq(maxval*1.1, (maxval*(1.0+(0.1*length(combtab)))), 0.1*maxval))
        if (length(combtab) > 0) { # failsafe for if there arent actually nay viable comps
            ## Failsafe for single comp and fix rounding error:
            failsafesignparam <- sign(maxval+(0.1*length(combtab)*valrange) - maxval*1.1)
            pout <- pout + stat_compare_means(comparisons = combtab, method = testtypeparam,
                                              # method.args = list(exact = TRUE), ### FOR SOME REASON - THIS IS KILLING MY RUN TIME, so omitting for now
                                              label = "p.format", size = 3,
                                              # label.y = seq(maxval*1.1, (maxval*(1.0+(0.1*length(combtab)))), 0.1*maxval))
                                              # label.y = seq(maxval*1.1, maxval*(0.1*length(combtab)*valrange), 0.1*valrange))
                                              label.y = seq(maxval*1.1, maxval+(0.1*length(combtab)*valrange), failsafesignparam*0.1*valrange))
            ## bug here where the rounding screws it up... really only when you have one comp though
            # label.y = seq(ceiling(maxval*1.1, 2), ceiling(maxval+(0.1*length(combtab)*valrange),2), round(0.1*valrange,2)))
        }
        
        # pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2)),
        #                     breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5))
        # pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(0.1*length(combtab)*valrange, digits = 2)),
        #                                   breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5))
        pout <- pout + scale_y_continuous(limits = c(round(minval-(abs(0.2*minval)), digits = 2), round(maxval+(0.1*length(combtab)*valrange+0.1*valrange), digits = 2)),
                                          # breaks = seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5))
                                          breaks = round(seq(round(minval-(abs(0.2*minval)), digits = 2), round(maxval*(1.1+(0.1*length(combtab))), digits = 2),length.out = 5),2))
        
        # expand_scale(mult = c(1.2, round(maxval*(1.1+(0.1*length(combtab))), digits = 2)))
        
    }
    
    ## Add in point coloring if desired # added 9/21/2021
    if (pointcoloring != TRUE) {
        pout <- pout + scale_color_manual(values = "black", breaks = "black")
    }
    
    pout
    
    #return(pout)
    
}




# --------------------------------- histogram plot ---------------------------------

## Histogram function
# p1datatab = as.data.frame(colSums(counttab))
# colnames(p1datatab) = "sampreadsum"
# plot_histogram(data = p1datatab,
#                fitcurve=TRUE, binwidthparam = NULL, limitx = NULL,
#                labsparam = list(title = "Counts for total reads per Sample prefilter", x = "Count", y = "reads in sample"),
#                pdffile = paste(outfilepath, "prefiltered_samp_readcount_hist.pdf", sep=""))
# plot_histogram <- function(data, fitcurve = FALSE, binwidthparam=NULL, labsparam, limitx = NULL, pdffile) {
plot_histogram <- function(data, groupvar = NULL, fitcurve = FALSE, binparam=NULL, labsparam, limitx = NULL, limity = NULL) {

    ## Amazingly - i need a catch here to make it a data.frame 2022-10-10
    data <- data.frame(data)
    
    ## If we want to fit curve to our single curve
    if (fitcurve == TRUE){
        meantemp = mean(data[,1])
        sdtemp = sd(data[,1])
    }

    ## binning data properly
    if (is.null(binparam)) {binparam = 10}
    # if (binparam > 100) {binparam = 100}  # 2022-03-10

    ## Add histograms for multiple groups
    if (!is.null(groupvar)) {
      ## Failsafe to make sure group is a character
      data$groupvar <- as.character(groupvar[,1])
    }
  
    ## plot the stuff
    pout <- ggplot(data, aes(x=data[,1]))
    if (!is.null(groupvar)) {
        pout <- pout + aes(fill=data$groupvar)
        pout <- pout + geom_histogram(color = "black", bins = binparam, alpha = 0.3, position="identity")
    } else {
        pout <- pout + geom_histogram(fill = "blue", color = "black", bins = binparam, alpha = 0.3, position="identity")
    }
    if (fitcurve == TRUE){
        pout <- pout + geom_freqpoly(mapping = aes(color = data$groupvar), bins = binparam, position = "identity")
        # pout <- pout + stat_function(fun = function(x, mean, sd) {
        #     dnorm(x=x, mean = meantemp, sd = sdtemp)*binparam}, color="red")
    }
    if (!is.null(limitx)) {pout <- pout + coord_cartesian(xlim = limitx)}
    if (!is.null(limity)) {pout <- pout + coord_cartesian(ylim = limity)}
    #pout <- pout + labs(title = "Counts for rowMeans per Gene prefilter", x = "Count", y = "log2(RowMean)")
    pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y)
    pout <- pout + theme_bw()
    

    return(pout)
}

# --------------------------------- PCA plotter ---------------------------------

## PCA plotting
# colorvar = metatablefilt3[,desccol,drop=FALSE]
# labsparam = c(list(title = "PCA", x = "PC1", y = "PC2", color=colnames(metatablefilt3)[desccol]))
# outfile = paste(outfilepath, "pca_plots/", colnames(metatablefilt3)[desccol], "_pca_plot.pdf", sep="")
# pca_plotter(pcadata = pcadata, colorvar = colorvar, scalecolor=FALSE, labelpoints = TRUE, labsparam = labsparam, outfile = outfile)

pca_plotter <- function(pcadata, colorvar=NULL, shapevar = NULL,
                        scalecolor=FALSE, labelpoints = FALSE, separatelegend = FALSE, labsparam,
                        returnoutliers = TRUE, custompcs = NULL) {
    pcdata = prcomp(pcadata, scale = TRUE)
    if (is.null(custompcs)){
      pc1label <- "PC1"
      pc2label <- "PC2"
    } else{
      pc1label <- custompcs[1]
      pc2label <- custompcs[2]
    }
    pc1 = pcdata$x[,pc1label]
    pc2 = pcdata$x[,pc2label]


    ## Add in outlier analysis
    pca_outliers <- NULL
    if (returnoutliers == TRUE) {
        zscore_pc <- scale(pc1)
        outliercutoff <- 3
        pca_outliers <- zscore_pc[zscore_pc > outliercutoff,,drop=FALSE]
        colnames(pca_outliers) <- "pca_diff"
    }

    pcloadings = as.data.frame(pcdata$rotation)
    pc1genes = rownames(pcloadings[order(pcloadings[,pc1label], decreasing=TRUE),])[1:min(ncol(pcadata), 5)]
    pc2genes = rownames(pcloadings[order(pcloadings[,pc2label], decreasing=TRUE),])[1:min(ncol(pcadata), 5)]
    pcproportions = summary(pcdata)$importance[2,]

    xlabel = paste(labsparam$x, "_", percent(pcproportions[pc1label]), "_", paste(pc1genes, collapse = "_"), sep="")
    ylabel = paste(labsparam$y, "_", percent(pcproportions[pc2label]), "_", paste(pc2genes, collapse = "_"), sep="")

    pcaplotdata = cbind.data.frame(PC1 = pc1, PC2 = pc2)
    if (!is.null(colorvar)) {pcaplotdata$colorvar = colorvar[,1]}
    if (!is.null(shapevar)) {pcaplotdata$shapevar = shapevar[,1]}

    if (labelpoints == TRUE) {
        pout <- ggplot(data = pcaplotdata, aes(x=pcaplotdata[,1], y = pcaplotdata[,2], label=rownames(pcaplotdata)))
    } else {
        pout <- ggplot(data = pcaplotdata, aes(x=pcaplotdata[,1], y = pcaplotdata[,2]))
    }
    # if (!is.null(colorvar)) {
    #   pout <- pout + aes(color=colorvar)
    #   pout <- pout + geom_point(size = 2.5, shape = 16)
    # }
    # if (!is.null(shapevar)) {
    #   pout <- pout + aes(shape=shapevar)
    #   pout <- pout + geom_point(size = 2.5)
    # }
    if (!is.null(colorvar)) {pout <- pout + aes(color=colorvar)}
    if (!is.null(shapevar)) {pout <- pout + aes(shape=shapevar)}
    if (!is.null(shapevar)) {
      pout <- pout + geom_point(alpha = 0.7)
    } else {
      pout <- pout + geom_point(alpha = 0.7, shape = 16)
    }
    
    if (labelpoints == TRUE) {pout <- pout + geom_text_repel(box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"), 
                                                             max.overlaps  = 1000)}
    if (returnoutliers == TRUE){
        pcaplotoutliers <- pcaplotdata[rownames(pca_outliers),]
        # pout <- pout + geom_text_repel(data = pcaplotoutliers, aes(x=pcaplotoutliers[,1], y = pcaplotoutliers[,2],
        #                                                            color=pcaplotoutliers[,3], label=rownames(pcaplotoutliers)),
        #                                                            box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
        if (!is.null(colorvar)) {
          pout <- pout + geom_text_repel(data = pcaplotoutliers, aes(x=pcaplotoutliers[,1], y = pcaplotoutliers[,2],
                                                                     color=pcaplotoutliers[,3], label=rownames(pcaplotoutliers)),
                                         box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
        } else {
          pout <- pout + geom_text_repel(data = pcaplotoutliers, aes(x=pcaplotoutliers[,1], y = pcaplotoutliers[,2], 
                                                                     label=rownames(pcaplotoutliers)),
                                         box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
        }
    }
    pout <- pout + labs(title = labsparam$title, x = xlabel, y = ylabel, color=labsparam$color, shape=labsparam$shape)
    if (scalecolor == TRUE) {pout <- pout + scale_color_gradient(low = "yellow", high = "blue")}

    if (separatelegend == TRUE) {
        legend <- g_legend(pout)
        pout <- pout + theme_pubr(legend = "none")
        return(list(pca_plot = pout, pca_legend = legend, pca_outliers = pca_outliers, pcdata = pcdata))
    } else {
        pout <- pout + theme_pubr(legend = "top")
        return(list(pca_out = pout, pca_outliers = pca_outliers, pcdata = pcdata))
    }

}

# --------------------------------- scatter plot ---------------------------------

## Scatter Plot
# input table has first two columns as x and y, any additional shapevar or colorvar must be a specified separate column
# datalabels work by giving a list of samples that match rownames that you want to label, NEED ROWNAMES IF YOU WANNA LABEL
# indata = data.frame(col1, col2, col3, col4)
# datalabels = rownames(intab[,col3 = LABELTHESE])
# colorvar = indata[,col3,drop=FALSE]
# shapevar = indata[,col4,drop=FALSE]
# labsparam = list(title = "raw to UQ norm count comparison", x = "raw", y = "norm")
# scatter_plotter(indata, colorvar, shapevar, datalabels, labsparam, plotstats)

scatter_plotter <- function(indata, colorvar=NULL, shapevar = NULL, sizevar = NULL, 
                            datalabels = NULL, labsparam, plotstats = FALSE, addjitterparam = NULL){

    if (!is.null(colorvar)) {indata$colorvar = colorvar[,1]}
    if (!is.null(shapevar)) {indata$shapevar = shapevar[,1]}
    if (!is.null(sizevar)) {indata$sizevar = sizevar[,1]}

    plotdata <- na.omit(indata)

    pout <- ggplot(data = plotdata, aes(x=plotdata[,1], y = plotdata[,2]))
    # if (!is.null(colorvar)) {
    #     pout <- pout + aes(color=plotdata$colorvar,fill=plotdata$colorvar)
    #     if (sum(areColors(plotdata$colorvar)) == length(plotdata$colorvar) & !is.numeric(plotdata$colorvar)) { ## If your colorvar are literally colors and NOT number - then use them
    #         pout <- pout + scale_color_identity()
    #     }
    # } ##1/13/2021
    if (!is.null(colorvar)) {
      pout <- pout + aes(color = colorvar,fill = colorvar)
      if (sum(areColors(plotdata$colorvar)) == length(plotdata$colorvar) & !is.numeric(plotdata$colorvar)) { ## If your colorvar are literally colors and NOT number - then use them
        pout <- pout + scale_color_identity()
        pout <- pout + scale_fill_identity()
      }
    } ##1/13/2021
    if (!is.null(shapevar)) {pout <- pout + aes(shape=plotdata$shapevar)}
    if (!is.null(sizevar)) {pout <- pout + aes(size=plotdata$sizevar)}
    # pout <- pout + geom_point(alpha = 0.7, shape = 16)
    if (!is.null(addjitterparam)){
      if (!is.null(shapevar)) {
        pout <- pout + geom_point(alpha = 0.7, position = position_jitter(width = addjitterparam))
      } else {
        pout <- pout + geom_point(alpha = 0.7, shape = 16, position = position_jitter(width = addjitterparam))
      }
    } else { ## Have to keep shapevar in their even if with jiutterparam
      if (!is.null(shapevar)) {
        pout <- pout + geom_point(alpha = 0.7)
      } else {
        pout <- pout + geom_point(alpha = 0.7, shape = 16)
      }
    }
    if (!is.null(datalabels)) {

        ## pout <- pout + geom_label_repel(data = labeledplottab, mapping = aes(x=labeledplottab[,1], y = labeledplottab[,2],
        ##                            # color=labeledplottab$colorvar, size=labeledplottab$sizevar,
        ##                            label=rownames(labeledplottab)),
        ##                            box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
        #pout <- pout + aes(label=rownames(labeledplottab))

        # labeledplottab <- plotdata[datalabels,]
        # pout <- pout + geom_label_repel(data = labeledplottab, mapping = aes(x=labeledplottab[,1], y = labeledplottab[,2],
        #                                 color = NULL, size = NULL, shape = NULL, label = rownames(labeledplottab)),
        #                                 box.padding = unit(1, "lines"), point.padding = unit(1, "lines"))

        plotdata$labelvar <- ""
        plotdata[datalabels,"labelvar"] <- datalabels
        pout <- pout + aes(label=plotdata$labelvar)
        pout <- pout + geom_label_repel(mapping = aes(size = NULL, label = plotdata$labelvar, fill = NULL),
                                        box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"),
                                        size = 2, max.overlaps = 10000)

    }
    pout


    if (plotstats != FALSE) {
        if (plotstats == TRUE) {corrmethod = "spearman"} else {corrmethod = plotstats}
        # pout <- pout + geom_smooth(method = lm, color="black")
        # pout <- pout + geom_text(label = lm_eqn(plotdata[,c(1:2)]), parse=TRUE,
        #                          mapping = aes(x = 0.8*max(plotdata[,1]), y = 0.9*max(plotdata[,2])), size = 8)
        #
        pout <- pout + geom_smooth(method = "gam", formula = y ~ x, color="black", data = plotdata[,c(1,2)], 
                                   inherit.aes = FALSE, aes(x = plotdata[,1], y = plotdata[,2], color = NULL, shape = NULL),)
        pout <- pout + stat_cor(data = plotdata[,c(1,2)], inherit.aes = FALSE, aes(x = plotdata[,1], y = plotdata[,2], color = NULL, shape = NULL),
                                method = corrmethod, label.x = 0.6*min(plotdata[,1]), label.y = 1.1*max(plotdata[,2]), size = 3)
    }
    pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, color = labsparam$color, shape = labsparam$shape, size = labsparam$size)
    pout <- pout + theme_pubr()
    pout
    
    return(pout)
}




# --------------------------------- scatter plot with BIG add on rasterizing ---------------------------------

# scatter_BIG_plotter <- function(indata, colorvar=NULL, shapevar = NULL, sizevar = NULL, 
#                             datalabels = NULL, labsparam, plotstats = FALSE, addjitterparam = NULL){
#     
#     if (!is.null(colorvar)) {indata$colorvar = colorvar[,1]}
#     if (!is.null(shapevar)) {indata$shapevar = shapevar[,1]}
#     if (!is.null(sizevar)) {indata$sizevar = sizevar[,1]}
#     
#     plotdata <- na.omit(indata)
#     
#     pout <- ggplot(data = plotdata, aes(x=plotdata[,1], y = plotdata[,2]))
#     # if (!is.null(colorvar)) {
#     #     pout <- pout + aes(color=plotdata$colorvar,fill=plotdata$colorvar)
#     #     if (sum(areColors(plotdata$colorvar)) == length(plotdata$colorvar) & !is.numeric(plotdata$colorvar)) { ## If your colorvar are literally colors and NOT number - then use them
#     #         pout <- pout + scale_color_identity()
#     #     }
#     # } ##1/13/2021
#     if (!is.null(colorvar)) {
#         pout <- pout + aes(color = colorvar,fill = colorvar)
#         if (sum(areColors(plotdata$colorvar)) == length(plotdata$colorvar) & !is.numeric(plotdata$colorvar)) { ## If your colorvar are literally colors and NOT number - then use them
#             pout <- pout + scale_color_identity()
#             pout <- pout + scale_fill_identity()
#         }
#     } ##1/13/2021
#     if (!is.null(shapevar)) {pout <- pout + aes(shape=plotdata$shapevar)}
#     if (!is.null(sizevar)) {pout <- pout + aes(size=plotdata$sizevar)}
#     # pout <- pout + geom_point(alpha = 0.7, shape = 16)
#     if (!is.null(addjitterparam)){
#         if (!is.null(shapevar)) {
#             pout <- pout + geom_point(alpha = 0.7, position = position_jitter(width = addjitterparam))
#         } else {
#             pout <- pout + geom_point(alpha = 0.7, shape = 16, position = position_jitter(width = addjitterparam))
#         }
#     } else { ## Have to keep shapevar in their even if with jiutterparam
#         if (!is.null(shapevar)) {
#             pout <- pout + geom_point(alpha = 0.7)
#         } else {
#             pout <- pout + geom_point(alpha = 0.7, shape = 16)
#         }
#     }
#     if (!is.null(datalabels)) {
#         
#         ## pout <- pout + geom_label_repel(data = labeledplottab, mapping = aes(x=labeledplottab[,1], y = labeledplottab[,2],
#         ##                            # color=labeledplottab$colorvar, size=labeledplottab$sizevar,
#         ##                            label=rownames(labeledplottab)),
#         ##                            box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
#         #pout <- pout + aes(label=rownames(labeledplottab))
#         
#         # labeledplottab <- plotdata[datalabels,]
#         # pout <- pout + geom_label_repel(data = labeledplottab, mapping = aes(x=labeledplottab[,1], y = labeledplottab[,2],
#         #                                 color = NULL, size = NULL, shape = NULL, label = rownames(labeledplottab)),
#         #                                 box.padding = unit(1, "lines"), point.padding = unit(1, "lines"))
#         
#         plotdata$labelvar <- ""
#         plotdata[datalabels,"labelvar"] <- datalabels
#         pout <- pout + aes(label=plotdata$labelvar)
#         pout <- pout + geom_label_repel(mapping = aes(size = NULL, label = plotdata$labelvar, fill = NULL),
#                                         box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"),
#                                         size = 2, max.overlaps = 10000)
#         
#     }
#     pout
#     
#     
#     if (plotstats != FALSE) {
#         if (plotstats == TRUE) {corrmethod = "spearman"} else {corrmethod = plotstats}
#         # pout <- pout + geom_smooth(method = lm, color="black")
#         # pout <- pout + geom_text(label = lm_eqn(plotdata[,c(1:2)]), parse=TRUE,
#         #                          mapping = aes(x = 0.8*max(plotdata[,1]), y = 0.9*max(plotdata[,2])), size = 8)
#         #
#         pout <- pout + geom_smooth(method = "gam", formula = y ~ x, color="black", data = plotdata[,c(1,2)], 
#                                    inherit.aes = FALSE, aes(x = plotdata[,1], y = plotdata[,2], color = NULL, shape = NULL),)
#         pout <- pout + stat_cor(data = plotdata[,c(1,2)], inherit.aes = FALSE, aes(x = plotdata[,1], y = plotdata[,2], color = NULL, shape = NULL),
#                                 method = corrmethod, label.x = 0.6*min(plotdata[,1]), label.y = 1.1*max(plotdata[,2]), size = 3)
#     }
#     pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, color = labsparam$color, shape = labsparam$shape, size = labsparam$size)
#     pout <- pout + theme_pubr()
#     pout
#     
#     return(pout)
# }


# --------------------------------- linechart plot ---------------------------------

## linechart Plot
# input table has first two columns as x and y, any additional shapevar or colorvar must be a specified separate column.
# datalabels work by giving a list of samples that match rownames that you want to label
# indata = data.frame(col1, col2, col3, col4)
# datalabels = rownames(intab[,col3 = LABELTHESE])
# colorvar = indata[,col3,drop=FALSE]
# shapevar = indata[,col4,drop=FALSE]
# labsparam = list(title = "raw to UQ norm count comparison", x = "raw", y = "norm")
# scatter_plotter(indata, colorvar, shapevar, datalabels, labsparam, plotstats)

linechart_plotter <- function(indata, groupvar = NULL, colorvar=NULL, shapevar = NULL,
                              datalabels = NULL, showpoints = TRUE, labsparam){

    plotdata <- indata

    if (!is.null(groupvar)) {plotdata$groupvar = groupvar[,1]}
    if (!is.null(colorvar)) {plotdata$colorvar = colorvar[,1]}
    if (!is.null(shapevar)) {plotdata$shapevar = shapevar[,1]}

    pout <- ggplot(data = plotdata, aes(x=plotdata[,1], y = plotdata[,2]))
    if (!is.null(colorvar)) {pout <- pout + aes(color=plotdata$colorvar)}
    if (!is.null(shapevar)) {pout <- pout + aes(shape=plotdata$shapevar)}
    if (!is.null(groupvar)) {pout <- pout + aes(group=plotdata$groupvar)}
    if (showpoints == TRUE) { pout <- pout + geom_point(alpha = 0.3, shape = 16)}
    pout <- pout + geom_path()
    if (!is.null(datalabels)) {
        labeledplottab <- plotdata[datalabels,]
        pout <- pout + geom_text_repel(data = labeledplottab, aes(x=labeledplottab[,1], y = labeledplottab[,2],
                                                                  color=labeledplottab$colorvar, label=rownames(labeledplottab)),
                                       box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
    }
    pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, color = labsparam$color, shape = labsparam$shape, size = labsparam$size)
    pout <- pout + theme_pubr()

    return(pout)
}






## Volcano Scatter Plot
# intable = table with log2fc as 2nd columns and pvalue as 1st
# pvalcutoff = value for pval cutoff
# log2fccutoff = value for log2fc cutoff
# labelparam = param for labeling plots, true or false - will label ALL GOI for now
# pout = create_volcano_plot(intable, pvalcutoff = 0.01, log2fccutoff = 1, labelparam = FALSE)
create_volcano_plot <- function(intable, pvalcutoff = 0.01, log2fccutoff = 1, labeledgenes = NULL, nameparam = "Volcano Plot") {
    plottab = intable
    plottab[is.na(plottab[,2]),2] <- 1 ## Failsafe to have any pvalues of NA be set as 1

    ## Failsafe for the log2fc - setting all NA to 0
    plottab[is.na(plottab[,1]),1] <- 0

    plottab$color = "notsig"
    plottab[plottab[,2] < pvalcutoff & plottab[,1,drop=FALSE] > 0, "color"] <- "sigpvalup" # pink
    plottab[plottab[,1,drop=FALSE] > log2fccutoff, "color"] <- "siglog2fcup" # red
    plottab[plottab[,1,drop=FALSE] > log2fccutoff & plottab[,2] < pvalcutoff, "color"] <- "sigbothup" # darkred

    plottab[plottab[,2] < pvalcutoff & plottab[,1,drop=FALSE] < 0, "color"] <- "sigpvaldown" # light blue
    plottab[plottab[,1,drop=FALSE] < -log2fccutoff, "color"] <- "siglog2fcdown" # blue
    plottab[plottab[,1,drop=FALSE] < -log2fccutoff & plottab[,2] < pvalcutoff, "color"] <- "sigbothdown" # dark blue

    plottab[,2] = -log10(plottab[,2])
    plottabin = as.data.frame(plottab)
    # great, pink, red, darkred, light blue, blue, dark blue
    colorparam = c("grey", "#ff9999", "#ff0000", "#990000", "#b2b2ff", "#1919ff", "#0000b2")
    names(colorparam) = c("notsig","sigpvalup", "siglog2fcup", "sigbothup", "sigpvaldown", "siglog2fcdown", "sigbothdown" )

    if (is.logical(labeledgenes) && labeledgenes == TRUE) {
        maxnumlabeledgenes = 40
        labeledgenes = rownames(plottabin[plottabin[,"color"] %in% c("sigbothup", "sigbothdown"),])
        labelmetric = plottabin[,1, drop=FALSE] * plottabin[,2, drop=FALSE]
        labeledgenes = rownames(labelmetric[order(abs(labelmetric), decreasing = TRUE),,drop=FALSE][1:40,,drop=FALSE])
        }
    labeltab = plottabin[rownames(plottabin) %in% labeledgenes,]

    pout <- ggplot(plottabin, aes(x = plottabin[,1], y = plottabin[,2], color = plottabin[,3]))
    pout <- pout + geom_point(size = 0.5, shape = 16)
    pout <- pout + scale_color_manual(values = colorparam)
    pout <- pout + labs(x = "log2fc", y = "-log10(pvalue)",
                        title = nameparam, color = "significance")
    if (!is.null(labeledgenes)) {pout <- pout + geom_label_repel(
        data = labeltab,
        aes(x = labeltab[,1], y=labeltab[,2], label = rownames(labeltab), color=labeltab[,3]),
        size=2, segment.size = 0.2, fontface="bold", segment.color = "black",
        box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"), 
        max.overlaps = 10000
    )}

    return(pout)

}




## WRITE HEATMAP FUNCTION
#counttab ## the counttable
#subsetnum = FALSE ## set to a number if you want to select by the N most varied genes (good for subsetting)
#metatable ## The metatable
#annotationlist ## The annotation object created beforehand (NEED TO AUTOMATE THIS MORE)
#colclusterparam = FALSE ## if not FALSE - will cluster the columns
#rowclusterparam = FALSE ## if not FALSE - will cluster the rows
#pdfoutfile ## the path to the outfile that the pdf will be saved to
# clustermethod = "ward.D2"
create_heatmap <- function(counttab,
                           subsetnum = FALSE,
                           scale_data = TRUE,
                           colmetatable = NULL, colannotationlist = NULL,
                           rowmetatable = NULL, rowannotationlist = NULL,
                           colclusterparam = FALSE, rowclusterparam = FALSE,
                           separate_legend = FALSE, heatmapcolorparam = NULL,
                           addborders = FALSE, columnsplitparam = NULL, rowsplitparam = NULL,
                           titleparam = c(column_title = "Samples", row_title = "Genes", figure_title = NULL, legend_title = "zscore"),
                           show_column_dend_param = TRUE, show_row_dend_param = TRUE ## 2023-01-30
                           # pdfoutfile
                           ) {

    colclusterparamIN = colclusterparam
    rowclusterparamIN = rowclusterparam

    ## Filter to the N most varied genes if applicable
    if (subsetnum !=FALSE) {
        counttab_variance = sort(apply(counttab,1,var), decreasing=TRUE)
        counttab = counttab[names(counttab_variance[1:subsetnum]),,drop=FALSE]
    }

    # ## Calc. spearman correlation and use values for column clustering before any other alterations
    #   if (colclusterparam != FALSE) {
    #   ## FAILSAFE - if there are any columns that sneak in somehow which have a variance of 0 - they have to be excluded for plotting
    #       nearZeroVarcols <- which(sapply(counttab, var) == 0)
    #       if (length(nearZeroVarcols) > 0) {
    #           counttab = counttab[,!colnames(counttab) %in% names(nearZeroVarcols)]
    #           print(paste("excluding column ", names(nearZeroVarcols), " due to zero variance", sep=""))}
    #       cordata <- cor(counttab, method="spearman")
    #       # coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
    #       coldistance = dist(t(as.matrix(na.omit(cordata))), method = "euclidean")
    #       colcluster = hclust(coldistance, method = clustermethod)
    #       colclusterparam = colcluster
    #   }
    
    ## Need failsafes for column and row clustering
    # if (is.logical(colclusterparamIN)) {
    colnearZeroVarcols <- rownearZeroVarcols <- NULL
    if (is.logical(colclusterparam)) {
        if (colclusterparam != FALSE) {
          colnearZeroVarcols <- which(apply(counttab, 2, var) == 0)
          if (length(colnearZeroVarcols) > 0) {
            counttab = counttab[,!colnames(counttab) %in% names(colnearZeroVarcols)]
            colmetatable <- colmetatable[!rownames(colmetatable) %in% names(colnearZeroVarcols),,drop=FALSE]
            print(paste("excluding column ", names(colnearZeroVarcols), " due to zero variance", sep=""))}
        }
    }
    # if (is.logical(rowclusterparamIN)) {
    if (is.logical(rowclusterparam)) {
        if (rowclusterparam != FALSE) {
            rownearZeroVarcols <- which(apply(counttab, 1, var) == 0)
            if (length(rownearZeroVarcols) > 0) {
                counttab = counttab[!rownames(counttab) %in% names(rownearZeroVarcols),]
                rowmetatable = rowmetatable[!rownames(rowmetatable) %in% names(rownearZeroVarcols),,drop=FALSE]
                print(paste("excluding row ", names(rownearZeroVarcols), " due to zero variance", sep=""))}
        }
    }
    
    ## Zscore out counttable, or turn into spearman correlation if doing sample-sample comparison
    # if (samplesample !=FALSE) {
    #     maptab <- cor(counttab, method="spearman")
    # } else {
        if(ncol(counttab) > 1 & scale_data == TRUE){
            maptab = as.matrix(t(apply(counttab, 1, function(x) zscore(x))))
        } else {
            maptab = as.matrix(counttab)
        }
        if (scale_data == "Log10") {
            maptab = as.matrix(log10(counttab+1))
        }
        if (scale_data == "Log2") {
            maptab = as.matrix(log2(counttab+1))
        }
    # }

    # if (rowclusterparam != FALSE) {
    #     rowdistance = dist(maptab, method = "euclidean")
    #     rowcluster = hclust(rowdistance, method = clustermethod)
    #     rowclusterparam = rowcluster
    # }

    ## Build Annotations from our metatable
    hatop = NULL
    if (!is.null(colannotationlist) & !is.null(colmetatable)) {
        temp1 <- vector("list", length(colannotationlist))
        names(temp1) = names(colannotationlist)
        annotlegendlist = lapply(temp1, function(x) x[[1]] =
                                     list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
        
        # ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
        for (annotnum in seq_len(length(annotlegendlist))) {
          annotname_select <- names(annotlegendlist)[annotnum]
          if (!is.null(names(colannotationlist[[annotname_select]]))) {
              annotlegendlist[[annotnum]] <- c(annotlegendlist[[annotnum]], list(at = names(colannotationlist[[annotname_select]])))
          } else {next}
        }
        ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
        
    ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
        showlegendparam = unname(unlist(lapply(colannotationlist, function(x) {
            numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
            is.null(numterms) || numterms <= 10})))
        hatop = HeatmapAnnotation(df = colmetatable,
                              col = colannotationlist,
                              na_col = "white",
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                              annotation_name_side = "left",
                              simple_anno_size = unit(min(60/length(colannotationlist), 5),"mm"),
                              show_legend = !separate_legend,
                              annotation_legend_param = annotlegendlist)
    }

    ## Defune the side annotation if data is supplied
    if (!is.null(rowannotationlist) & !is.null(rowmetatable)) {

        ## Define parameters for each of the labels on the annotation bars
        temp1 <- vector("list", length(rowannotationlist))
        names(temp1) = names(rowannotationlist)
        annotlegendlist = lapply(temp1, function(x)
            x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))
        
        ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order
        for (annotnum in seq_len(length(annotlegendlist))) {
            annotname_select <- names(annotlegendlist)[annotnum]
            if (!is.null(names(colannotationlist[[annotname_select]]))) {
                annotlegendlist[[annotnum]] <- c(annotlegendlist[[annotnum]], list(at = names(rowannotationlist[[annotname_select]])))
            } else {next}
        }
        ## ADDED 03-24-2022 - keeps the order of the annotationlist for the legend order

        ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
        showlegendparam = unname(unlist(lapply(rowannotationlist, function(x) {
            numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
            is.null(numterms) || numterms <= 10})))

        ## Look for any empty annotations - fill them with white, and later, make sure to hide their legend
        emptyannots = names(sapply(rowannotationlist, length)[sapply(rowannotationlist, length)==0])
        if (length(emptyannots) > 0){
            for (i in 1:length(emptyannots)) {
                temp1 = "white"
                names(temp1) = emptyannots[i]
                rowannotationlist[[emptyannots[i]]] = temp1
            }
            showlegendparam[which(names(rowannotationlist) %in% emptyannots)] = FALSE
        }
        ## Add param that will bolden the side annotation bars if it is <100, and omit the grid lines if more
        if (nrow(rowmetatable) < 100) {
            sideannotation_linebold_param <- gpar(fontsize = 0.5)} else {sideannotation_linebold_param <- NULL}
        haside = rowAnnotation(df = rowmetatable,
                           col = rowannotationlist,
                           na_col = "white",
                           # gp = gpar(fontsize = 0.01),
                           gp = sideannotation_linebold_param,
                           show_annotation_name=TRUE,
                           annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                           annotation_name_side = "top",
                           simple_anno_size = unit(min(60/length(rowannotationlist), 5),"mm"),
                           show_legend = !separate_legend,
                           annotation_legend_param = annotlegendlist)
    }

    ## Define the Heatmap
    # if (scale_data == TRUE) {
    #     if (samplesample == FALSE) {
    #         heatmapcolorparam = colorRamp2(c(-3,0,3), c("blue", "white", "red"))
    #     } else {
    #         heatmapcolorparam = colorRamp2(c((min(maptab)),1), c("white", "red"))
    #     }
    # } else {
    #     heatmapcolorparam = NULL
    # }

    lowcolor = midcolor = highcolor = lowvalue = midvalue = highvalue = NULL
    if(max(maptab, na.rm = TRUE) > 0) {
        highcolor = "red"
        highvalue = max(max(maptab, na.rm = TRUE), 0.2)
    } else {
        highcolor = "white"
        highvalue = 0
    }
    if(min(maptab, na.rm = TRUE) < 0) {
        lowcolor = "blue"
        lowvalue = min(min(maptab, na.rm = TRUE), -0.2)
    } else {
        lowcolor = "white"
        lowvalue = 0
    }
    if(max(maptab, na.rm = TRUE) > 0 & min(maptab, na.rm = TRUE) < 0) {
        midcolor = "white"
        midvalue = 0
    }
    if(is.null(heatmapcolorparam)) {
        heatmapcolorparam <- colorRamp2(
            c(lowvalue, midvalue, highvalue), c(lowcolor, midcolor, highcolor))
    } else{
        heatmapcolorparam <- heatmapcolorparam
    }
    
    # By default, dendrogram reordering is turned on if cluster_rows/cluster_columns is set as logical value or a clustering function. It is turned off if cluster_rows/cluster_columns is set as clustering object.
    
    ## This is so wasteful, but the only way to get borders to be toggleable is to repeat the code, ugh, oh well
    if (addborders == FALSE) {
      ht1 = Heatmap(maptab,
                    col = heatmapcolorparam, row_title = unname(titleparam["row_title"]), column_title = unname(titleparam["column_title"]), # Add in param to adjust titles 2023-01-25
                    cluster_columns = colclusterparamIN, cluster_rows = rowclusterparamIN,
                    clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                    #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
                    # 2022-08-25
                    column_dend_reorder = TRUE,
                    row_dend_reorder = TRUE,
                    # 2022-08-25
                    ## Add in columnsplitparam
                    column_split = columnsplitparam[colnames(maptab), ], # Add in care samples or features are removed due to no var - 2023-01-20
                    cluster_column_slices = FALSE,
                    row_split = rowsplitparam[rownames(maptab), ], # Add in care samples or features are removed due to no var - 2023-01-20
                    
                    ## Add borders
                    border = FALSE,
                    
                    show_column_names = ncol(maptab) <=100, column_names_gp = gpar(fontsize = 6),
                    show_row_names = nrow(maptab) <=200, row_names_side = "left", row_names_gp = gpar(fontsize=6),
                    show_row_dend = show_row_dend_param, show_column_dend = show_column_dend_param,
                    
                    heatmap_legend_param = list(
                      # title = ifelse(samplesample==FALSE, "Zscore", "Spearman\nCorrelation"),
                      title = titleparam["legend_title"],
                      # legend_height = unit(2.5, "cm"),
                      title_gp = gpar(fontsize = 8, fontface = "bold")),
                    top_annotation = hatop,
                    height = unit(min((nrow(maptab)/2), 12),"cm"), width = unit(min(ncol(maptab), 18),"cm")
                    
      )
    } else {
      ht1 = Heatmap(maptab,
                    col = heatmapcolorparam, row_title = unname(titleparam["row_title"]), column_title = unname(titleparam["column_title"]), # Add in param to adjust titles 2023-01-25
                    cluster_columns = colclusterparamIN, cluster_rows = rowclusterparamIN,
                    clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                    #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
                    
                    ## Add in columnsplitparam
                    column_split = columnsplitparam[colnames(maptab), ], # Add in care samples or features are removed due to no var - 2023-01-20
                    cluster_column_slices = FALSE,
                    row_split = rowsplitparam[rownames(maptab), ], # Add in care samples or features are removed due to no var - 2023-01-20
                    
                    ## Add borders
                    border = TRUE,
                    rect_gp = gpar(col = "black", lwd = 0.5),
                    
                    show_column_names = ncol(maptab) <=100, column_names_gp = gpar(fontsize = 6),
                    show_row_names = nrow(maptab) <=200, row_names_side = "left", row_names_gp = gpar(fontsize=6),
                    show_row_dend = show_row_dend_param, show_column_dend = show_column_dend_param,
                    
                    heatmap_legend_param = list(
                      # title = ifelse(samplesample==FALSE, "Zscore", "Spearman\nCorrelation"),
                      title = titleparam["legend_title"],
                      # legend_height = unit(2.5, "cm"),
                      title_gp = gpar(fontsize = 8, fontface = "bold")),
                    top_annotation = hatop,
                    height = unit(min((nrow(maptab)/2), 12),"cm"), width = unit(min(ncol(maptab), 18),"cm")
                    
      )
    }

    if (separate_legend == TRUE) {
        lgdoutlist = list()
        lgdinlist <- c(colannotationlist, rowannotationlist)
        # Create legendlist for col annots
        for (legendnum in seq_len(length(lgdinlist))){
            temp = lgdinlist[[legendnum]]
            if (is.function(temp)){
                lgdoutlist[[legendnum]] = Legend(col_fun = temp,
                                                 title = names(lgdinlist[legendnum]), legend_height = unit(2.5, "cm"),
                                                 title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
            } else {
                lgdoutlist[[legendnum]] = Legend(labels = names(temp), legend_gp = gpar(fill = unname(temp)),
                                                 title = names(lgdinlist[legendnum]), grid_height = unit(0.8, "cm"),
                                                 title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
            }
        }
        # for (legendnum in seq_len(length(colannotationlist))){
        #     temp = colannotationlist[[legendnum]]
        #     if (is.function(temp)){
        #         lgdoutlist[[legendnum]] = Legend(col_fun = temp,
        #                                        title = names(colannotationlist[legendnum]), legend_height = unit(2.5, "cm"),
        #                                        title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
        #     } else {
        #         lgdoutlist[[legendnum]] = Legend(labels = names(temp), legend_gp = gpar(fill = unname(temp)),
        #                                        title = names(colannotationlist[legendnum]), grid_height = unit(0.8, "cm"),
        #                                        title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
        #     }
        # }
        # # Create legendlist for row annots
        # for (legendnum in seq_len(length(rowannotationlist))){
        #     temp = rowannotationlist[[legendnum]]
        #     if (is.function(temp)){
        #         lgdoutlist[[legendnum]] = Legend(col_fun = temp,
        #                                          title = names(rowannotationlist[legendnum]), legend_height = unit(2.5, "cm"),
        #                                          title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
        #     } else {
        #         lgdoutlist[[legendnum]] = Legend(labels = names(temp), legend_gp = gpar(fill = unname(temp)),
        #                                          title = names(rowannotationlist[legendnum]), grid_height = unit(0.8, "cm"),
        #                                          title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
        #     }
        # }
        if (length(lgdoutlist) > 0) {
            heatmapannotationlegend = packLegend(
                list = lgdoutlist, max_height = unit(8.5, "in"), column_gap = unit(0.5, "cm"))
        } else {
            heatmapannotationlegend = NULL
        }
    } else {heatmapannotationlegend = NULL}

    #Legend()
    #annotlegendlist = lapply(temp1, function(x) x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))

    ## Plot out the heatmap
    # pdf(file = pdfoutfile, width=11,height=8.5)
    if (!is.null(rowannotationlist) & !is.null(rowmetatable)) {
        out1 <- draw(ht1 + haside, annotation_legend_side = "bottom", padding = unit(c(5, 20, 5, 5), "mm"), column_title = unname(titleparam["figure_title"]))
    } else {
        out1 <- draw(ht1, annotation_legend_side = "bottom", padding = unit(c(5, 20, 5, 5), "mm"), column_title = unname(titleparam["figure_title"]))
    }
    # grid.newpage()
    # if (!is.null(heatmapannotationlegend)) {draw(heatmapannotationlegend)}
    # junk <- dev.off()
    return(Filter(Negate(is.null), list(heatmap = out1, heatmaplegend = heatmapannotationlegend)))

}




## WRITE HEATMAP FUNCTION
#counttab ## the counttable
#subsetnum = FALSE ## set to a number if you want to select by the N most varied genes (good for subsetting)
#metatable ## The metatable
#annotationlist ## The annotation object created beforehand (NEED TO AUTOMATE THIS MORE)
#colclusterparam = FALSE ## if not FALSE - will cluster the columns
#rowclusterparam = FALSE ## if not FALSE - will cluster the rows
#pdfoutfile ## the path to the outfile that the pdf will be saved to
create_SS_heatmap <- function(counttab, metatable=NULL, annotationlist=NULL,
                              rowclusterparam=FALSE, colclusterparam=FALSE,
                              separate_legend = FALSE) {

    cordata <- cor(counttab, method="spearman")

    if (rowclusterparam != FALSE) {
        rowdistance = dist(t(as.matrix(cordata)), method = "euclidean")
        rowcluster = hclust(rowdistance, method = "ward.D2")
        rowclusterparam = rowcluster
    }
    if (colclusterparam != FALSE) {
        coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
        colcluster = hclust(coldistance, method = "ward.D2")
        colclusterparam = colcluster
    }

    ## Build Annotations from our metatable
    # hatop = HeatmapAnnotation(df = metatable,
    #                         col = annotationlist,
    #                         show_annotation_name = TRUE,
    #                         annotation_name_side = "left")

    ## Build Annotations from our metatable
    hatop = NULL
    if (!is.null(annotationlist) & !is.null(metatable)) {
        temp1 <- vector("list", length(annotationlist))
        names(temp1) = names(annotationlist)
        annotlegendlist = lapply(temp1, function(x) x[[1]] =
                                     list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
        ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
        showlegendparam = unname(unlist(lapply(annotationlist, function(x) {
            numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
            is.null(numterms) || numterms <= 10})))
        hatop = HeatmapAnnotation(df = metatable,
                                  col = annotationlist,
                                  na_col = "white",
                                  show_annotation_name = TRUE,
                                  annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                                  annotation_name_side = "left",
                                  simple_anno_size = unit(min(60/length(annotationlist), 5),"mm"),
                                  show_legend = !separate_legend,
                                  annotation_legend_param = annotlegendlist)
    }

    ## Define the Heatmap
    ht1 = Heatmap(t(as.matrix(cordata)),
                #col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),    ## Define the color scale for the heatmap
                row_title = "Samples",                                       ## Name the rows
                column_title = "Samples",                                  ## Name the columns

                cluster_columns = colclusterparam,                         ## Cluster the columns or leave as is
                cluster_rows = rowclusterparam,                            ## Cluster the rows or leave as is

                show_column_names = nrow(cordata) < 100,                                  ## Show the Column Names
                column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
                show_row_names = TRUE,                                    ## Show the row names
                row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
                row_names_gp = gpar(fontsize=6),

                show_row_dend = TRUE,                                     ## Show the dendrogram on the rows
                show_column_dend = TRUE,                                   ## Show the dendrogram on the columns

                heatmap_legend_param = list(title = "Spearman\nCorrelation",
                                            legend_height = unit(2, "cm"),
                                            title_gp = gpar(fontsize = 8)),
                top_annotation = hatop
    )

    ## Plot out the heatmap
    # pdf(file = pdfoutfile, width=11.5,height=8)
    # draw(ht1)
    # junk <- dev.off()

    if (separate_legend == TRUE) {
        lgdoutlist = list()
        for (legendnum in seq_len(length(annotationlist))){
            temp = annotationlist[[legendnum]]
            if (is.function(temp)){
                lgdoutlist[[legendnum]] = Legend(col_fun = temp,
                                                 title = names(annotationlist[legendnum]), legend_height = unit(2.5, "cm"),
                                                 title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
            } else {
                lgdoutlist[[legendnum]] = Legend(labels = names(temp), legend_gp = gpar(fill = unname(temp)),
                                                 title = names(annotationlist[legendnum]), grid_height = unit(0.8, "cm"),
                                                 title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))
            }
        }
        if (length(lgdoutlist) > 0) {
            heatmapannotationlegend = packLegend(
                list = lgdoutlist, max_height = unit(8.5, "in"), column_gap = unit(0.5, "cm"))
        } else {
            heatmapannotationlegend = NULL
        }
    } else {heatmapannotationlegend = NULL}

    #Legend()
    #annotlegendlist = lapply(temp1, function(x) x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))

    ## Plot out the heatmap
    # pdf(file = pdfoutfile, width=11,height=8.5)
    out1 <- draw(ht1, annotation_legend_side = "bottom", padding = unit(c(5, 20, 5, 5), "mm"))
    # grid.newpage()
    # if (!is.null(heatmapannotationlegend)) {draw(heatmapannotationlegend)}
    # junk <- dev.off()
    return(Filter(Negate(is.null), list(heatmap = out1, heatmaplegend = heatmapannotationlegend)))

}

## Input the metatable and this will build an annotation, if you want to preset
## the colors - then input them in proper list form, either as named list of
## named colors, or with color brewer
#' Create the annotation object for plotting in a heatmap
#' @param metatable the metatable containing information for the columns
#' @param customcolorlist DEFAULT: NULL, enter colorlist to manually set colors
#' @return return the annotation object
#' @keywords outliers
#' @import stats circlize RColorBrewer viridis
#' @export
#' @examples
#' metatable = data.frame(row.names = c("samp1", "samp2", "samp3", "samp4"),
#'     A = c(rep("high", 2), rep("low", 2)), B = seq(1,7,2))
#' customcolorlist = list(A = c("high" = "red", "low" = "blue"),
#'                        B = circlize::colorRamp2(seq(-5, 5, length = 3),
#'                        RColorBrewer::brewer.pal(3, "Reds")))
#' annotationlist_builder(metatable, customcolorlist)
annotationlist_builder <- function(metatable, customcolorlist = NULL){
  annotlist <- list()
  #colorlist = colors()[grep('gr(a|e)y', colors(), invert = TRUE)]
  ## Create the color list from manual distinct colors, expanded appropriately
  colorlist = rep(c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00",
    "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
    "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"),10)

  colorcount = 0
  for (colnum in 1:ncol(metatable)) {

    if (!is.null(customcolorlist) && colnames(metatable[,colnum,drop=FALSE])
        %in% names(customcolorlist)) {
      annotlist[colnames(metatable[,colnum,drop=FALSE])] =
        customcolorlist[colnames(metatable[,colnum,drop=FALSE])]}

    if (!is.null(customcolorlist) && !colnames(metatable[,colnum,drop=FALSE])
        %in% names(customcolorlist) | is.null(customcolorlist)) {

        ## Added failsafe for blank annotation columns
      if (sum(is.na(metatable[,colnum,drop=FALSE])) == nrow(metatable[,colnum,drop=FALSE])) {
          annotlist[[colnum]] = NA
          names(annotlist)[[colnum]] = colnames(metatable)[colnum]
      } else {

          if (!is.numeric(metatable[,colnum])) {

            #annotlist[[colnum]] = sample(colorlist, length(na.omit(unique(metatable[,colnum]))))
            annotlist[[colnum]] = colorlist[(1+colorcount):(colorcount + length(na.omit(unique(metatable[,colnum]))))]
            colorcount = colorcount + length(na.omit(unique(metatable[,colnum])))

            names(annotlist[[colnum]]) = na.omit(unique(metatable[,colnum]))
            names(annotlist)[[colnum]] = colnames(metatable)[colnum]
          }
          else {
              ## Add fail safe for numeric ranges of only one value
              if (min(metatable[,colnum], na.rm = TRUE) != max(metatable[,colnum], na.rm = TRUE)) {
              collist = colorRamp2(seq(min(metatable[,colnum], na.rm = TRUE),
                  max(metatable[,colnum], na.rm = TRUE), length = 3),
                  brewer.pal(3, "Purples"))
              } else {
                  singleval = min(metatable[,colnum], na.rm = TRUE)
                  if (singleval > 0) {
                      minval = 0
                      maxval = singleval
                  }
                  if (singleval < 0) {
                      minval = singleval
                      maxval = 0
                  }
                  if (singleval == 0) {
                      minval = singleval
                      maxval = 1
                  }
                  collist = colorRamp2(seq(minval, maxval, length = 3),
                                       brewer.pal(3, "Purples"))
              }
              annotlist[[colnum]] = collist
              names(annotlist)[colnum] = colnames(metatable)[colnum]
          }
      }
    }
  }

  return(annotlist)
}


# https://cran.r-project.org/web/packages/umap/vignettes/umap.html
## TSNE PLOT

# indatafile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run1c_rmoutliers2/rna_processing/normcounttab.txt"
# indata <- read.table(indatafile, sep = "\t", header = TRUE, row.names = 1)
# tsne_out <- Rtsne(t(indata))

## Well have clustering options as an input as well - which is a list of types of clustering, and then Ns of clustering
# clusteringparam = list(hierarchical = c(2,3,4,5), kmeans = c(2,3,4,5))
tsne_plotter <- function(indata, colorvar, labsparam, clusteringparam = NULL, seedparam = NULL, outfile) {
    
    ## Generate tsne data based off of the input data
    set.seed(ifelse(is.null(seedparam), 42, seedparam))
    # note that it has to be flipped to have samples on y-axis - THE INPUT OBJECT SHOULD ALREADY BE FLIPPED
    tsne_out <- Rtsne(indata)
    
    ## Create the plot object
    tsne_plotobject <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2])
    
    if (!is.null(clusteringparam)) {
        ## First add all of the columns I need
        tsne_plotobject <- cbind.data.frame(tsne_plotobject, t(unlist(clusteringparam)))
        colnames(tsne_plotobject)[3:(2+length(unlist(clusteringparam)))] <- 
            paste0(c(rep("hierarchical_cut_", length(clusteringparam$hierarchical)), 
                     rep("kmeans_cut_", length(clusteringparam$kmeans))), unlist(clusteringparam))
        
        ## Then for each type of clustering - add those as coloring params
        if ("hierarchical" %in% names(clusteringparam)) {
            fit_cluster_hierarchical=hclust(dist(tsne_plotobject[,c("x","y")]))
            for (hcutnum in seq_len(length(clusteringparam$hierarchical))) {
                selcol <- 2 + hcutnum
                tsne_plotobject[,selcol] <- factor(cutree(fit_cluster_hierarchical, k=clusteringparam$hierarchical[hcutnum]))
            }
        }
        if ("kmeans" %in% names(clusteringparam)) {
            
            for (kcutnum in seq_len(length(clusteringparam$hierarchical))) {
                selcol <- 2 + hcutnum + kcutnum
                fit_cluster_kmeans=kmeans(tsne_plotobject[,c("x","y")], clusteringparam$kmeans[kcutnum])
                tsne_plotobject[,selcol] <- factor(fit_cluster_kmeans$cluster)
            }
        }
    }
    
    if (!is.null(colorvar)) {tsne_plotobject$manualcolorvar = colorvar[,1]}
    # No way to keep the labels - but it does preserve the order, so just have to add back the labels
    rownames(tsne_plotobject) <- rownames(indata)
    
    ## Create the labels
    xlabel <- ifelse(is.null(labsparam$x), "TSNE_dimension_1", labsparam$x)
    ylabel <- ifelse(is.null(labsparam$y), "TSNE_dimension_2", labsparam$y)
    
    ### NEED TO RUN A FORLOOP ACROSS EACH OF THESE COLOR LABELS
    # Start the ggplot object
    
    numcolorparam <- ncol(tsne_plotobject) - 2
    plotoutlist <- list()
    varlist <- as.list(colnames(tsne_plotobject)[3:ncol(tsne_plotobject)])
    for (colorcolumnnum in seq_len(numcolorparam)) {
        colorcolsel <- colnames(tsne_plotobject)[2+colorcolumnnum]
        pout <- ggplot(data = tsne_plotobject, aes(x=tsne_plotobject[,1], y = tsne_plotobject[,2]))
        pout <- pout + aes(color=tsne_plotobject[,varlist[[colorcolumnnum]][1]])
        pout <- pout + geom_point(size = 2.5, shape = 16)
        pout <- pout + labs(title = labsparam$title, x = xlabel, y = ylabel, color=colorcolsel)
        pout <- pout + theme_pubr(legend = "top")
        plotoutlist[[colorcolumnnum]] = pout
        names(plotoutlist)[colorcolumnnum] <- colorcolsel
    }
    
    pdf(outfile, useDingbats = FALSE)
    for (colorcolumnnum in seq_len(numcolorparam)) {
        print(plotoutlist[[colorcolumnnum]])
    }
    junk <- dev.off()
    
    return(list(tsne_plot = pout, tsne_table = tsne_plotobject))
    
}






## Get legend function to print on a separate page
# g2 <- ggplot()
# g2 <- g2 + geom_bar(aes(y = value, x = sample, fill = gene), data = propplottab, stat="identity", size=0.5, color="black")
# g2 <- g2 + scale_fill_manual(values = plotcolors)
# legend <- g_legend(g2)
# g2 <- g2 + theme(legend.position = "none")
#
# pdf(paste(outfilepath, "read_distribution_bar_chart.pdf", sep=""))
# print(g2)
# grid.newpage()
# grid.draw(legend)
# dev.off()
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


## Cbind fill function
cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


## Zscore Function
zscore <- function(x) {
  y=(x-mean(x, na.rm=TRUE))/sd(x, na.rm = TRUE)
  return(y)
}

## lm_eqn function for stat plotting on scatterplots
lm_eqn = function(df){
    colnames(df) = c("x", "y")
    m = lm(y ~ x, df);
    # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
    #                  list(a = format(coef(m)[1], digits = 2),
    #                       b = format(coef(m)[2], digits = 2),
    #                       r2 = format(summary(m)$r.squared, digits = 3)))
    eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
    as.character(as.expression(eq));
}


## parallel processing template
# ## Due to the number of iterations that this code has to run through - it has the option to run in parallel, input the number of cores here. This can be checked using < detectCores() > after loading the library("parallel")
# cores = 6
# ## Start the parallel processing cluster with the given number of cores
# suppressMessages(library(doParallel))
# cl <- makeCluster(cores)
# registerDoParallel(cl)

# ## Define combining function for running in parallel
# comb <- function(x, ...) {lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))}
# boxplotfuncout <- foreach(dbcount=1:nrow(dbtab), .packages=packagelist, .combine = 'comb', .multicombine = TRUE, .init=list(list(), list(), list(), list())) %dopar% {
#   DO WORK
# }
# stopCluster(cl)









## Supp function to check if strings are colors:
## https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
areColors <- function(x) {
    sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), 
                 error = function(e) FALSE)
    })
}



## Cohens d for effect size
cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    # md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    md  <- abs(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))        ## edit added 9/20/2021
    # csd <- lx * var(x) + ly * var(y)
    csd <- lx * var(x, na.rm = TRUE) + ly * var(y, na.rm = TRUE)     ## edit added 9/20/2021
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation
    cd  <- md/csd                        ## cohen's d
    return(cd)
}


## Create a custom ggally function with coloring the correlation values
topcorr_fn <- function(data, mapping, method = "spearman", use="pairwise.complete.obs", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

custom_corrplot <- function(intable, methodparam = "spearman") {
  corplotout <- ggpairs(intable, 
                upper = list(continuous = topcorr_fn),
                lower = list(continuous = "smooth"))
  corplotout <- corplotout + theme()
  return(corplotout)
  
}



create_waffle_plot <- function(waffledata, colorparam = NULL, rowparam = NULL) {
  # Turn bartab1 into data for a waffle part
  # waffle_data <- unlist(apply(waffledata, 1, function(x) rep(x[1], x[2])))
  waffle_plottab <- waffledata[,2]
  names(waffle_plottab) <- waffledata[,1]
  
  # Set colors for the waffleplot
  colorlist <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00",
                 "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
                 "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                 "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown")
  if (is.null(colorparam)) {
    wafflecolors <- colorlist[1:length(waffle_plottab)]
  } else {
    wafflecolors <- colorparam
  }
  
  
  if (is.null(rowparam)) {rowparam <- ceiling(sum(waffle_plottab)/25)}
  waffle_out <- waffle(waffle_plottab, rows = rowparam, colors = wafflecolors, size = 0.5, xlab = "1 Square = 1 Student")
  waffle_out <- waffle_out + ggtitle(label = colnames(waffledata)[1])
  
  return(waffle_out)
  
}




# simpleCap <- function(x) {
#   s <- strsplit(x, " ")[[1]]
#   paste(toupper(substring(s, 1,1)), substring(s, 2),
#         sep="", collapse=" ")
# }

simpleCap <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}





# --------------------------------- boxplots for a counttable based off of annotation ---------------------------------

boxplots_from_counttable_by_annotation <- function(counttable, boxplot_annotationtable, outfilepath, 
                                                   calculate_corrvalues = FALSE, 
                                                   bp_color_ref = NULL, violinparam = FALSE, plotstatsparam = "inter", plotpointsparam = TRUE) {
    # These are the features and annotations we will run over
    bp_features <- rownames(counttable)
    bp_annots <- colnames(boxplot_annotationtable)
    
    params_grid <- expand.grid(bp_features, bp_annots)
    colnames(params_grid) <- c("bp_features", "bp_annots")
    wilcox_statoutlist <- spearman_statoutlist <- list()
    for (params_num in seq_len(nrow(params_grid))) {
        
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
        bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "feature",
                                 labsparam = list(title = paste0(feature_sel, " values for samples by ", annot_sel), x = annot_sel, y = feature_sel,
                                                  catorder = "cohort", featorder = featorder_param),
                                 plotstats = plotstatsparam, testtypeparam = "wilcox.test", violin = violinparam, plotpoints = plotpointsparam,
        )
        
        # Add custom coloring for this particular section
        # bp_color_ref <- custom_annotation_list_from_colorguide(COI=group_column_selected, colorguide)
        if (!is.null(bp_color_ref)){
            if (annot_sel %in% names(bp_color_ref)) {
                bpout <- bpout + scale_fill_manual(breaks = names(bp_color_ref[[annot_sel]]), values = unname(bp_color_ref[[annot_sel]]))
            }
        }
        # Add custom coloring for this particular section
        if (!violinparam) {
            pdf(paste0(outfilepath, annot_sel, "/", feature_sel, "_boxplot.pdf"), useDingbats = FALSE)
        } else {
            pdf(paste0(outfilepath, annot_sel, "/", feature_sel, "_violin.pdf"), useDingbats = FALSE)
        }
        
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


# --------------------------------- convenience function to grab limits from one plot to apply to another ---------------------------------

get_plot_limits <- function(plot) {
    gb = ggplot_build(plot)
    xmin = gb$layout$panel_params[[1]]$x.range[1]
    xmax = gb$layout$panel_params[[1]]$x.range[2]
    ymin = gb$layout$panel_params[[1]]$y.range[1]
    ymax = gb$layout$panel_params[[1]]$y.range[2]
    list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}
# get_plot_limits(pout1)
# save_plot_limits <- unlist(get_plot_limits(pout1))
# pout2 <- pout2 + coord_cartesian(xlim = save_plot_limits[c("xmin", "xmax")], ylim = save_plot_limits[c("ymin", "ymax")])



