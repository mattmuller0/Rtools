################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## WGCNA Analysis

## Load in libraries
packagelist = c("WGCNA", "tools", "blacksheepr", "devtools")
junk <- lapply(packagelist, function(xxx) suppressMessages(
    require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/mgc_plotting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/mgc_geneset_analysis_functions.R')

# # INPUT
# # counttable - the normalized counttable to derive modules from
# # metatable - the metatable to explore correlations against - will binarize any categorical variables
# ## Read in our platelet data
# incounttablefile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/covid-platelet-sle-hcq/output/run4-updatedinfo/rna_processing/normcounttab.txt"
# innormcounttab <- read.table(incounttablefile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
# 
# ## Read in metadata
# inmetafile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/covid-platelet-sle-hcq/output/run4-updatedinfo/rna_processing/metatable_in.csv"
# inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
# 
# # Subset our sle data for just our lupus pts
# slecounttab <- innormcounttab[,rownames(inmetatable[inmetatable[,"cohort"] == "sle",])]
# 
# # The intable needs to be in ML form with features as cols
# counttable <- as.data.frame(t(slecounttab))
# # metatable <- inmetatable[inmetatable[,"cohort"] == "sle",2:ncol(metatable)]
# metatable <- inmetatable[inmetatable[,"cohort"] == "sle",
#                          c("hcq", "on_hcq_type_drug", "dose..mg.", colnames(inmetatable)[grepl("comp_", colnames(inmetatable))])]
# # temp_meta <- metatable[,]
# 
# outfilepath <- "/Users/tosh/Desktop/Ruggles_Lab/projects/covid-platelet-sle-hcq/output/run4-updatedinfo/WGCNA/sle_only/"
# dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)

# # PREPROCESSING - probably will never need
# # first looking for genes and samples with too many missing values (should be none?)
# qccheck1 <- goodSamplesGenes(counttable)
# if (!gsg$allOK) {
#     # Optionally, print the gene and sample names that were removed:
#     if (sum(!gsg$goodGenes)>0)
#         printFlush(paste("Removing genes:", paste(names(counttable)[!gsg$goodGenes], collapse = ", ")));
#     if (sum(!gsg$goodSamples)>0)
#         printFlush(paste("Removing samples:", paste(rownames(counttable)[!gsg$goodSamples], collapse = ", ")));
#     # Remove the offending genes and samples from the data:
#     datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
# }

## Determine power and create network

create_WGCNA_modules <- function(counttable, powerparam = NULL, minModuleSizeParam = 30,
                                 savedendroplotfile = NULL) {
    
    ## Create outfolder in case it wasnt
    # dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)
    
    powers = c(1:20)
    sft = pickSoftThreshold(counttable, powerVector = powers, verbose = 5)
    ## Determine the power needed
    if (is.null(powerparam)) {
        powerthresh = sft$powerEstimate
        if (is.na(powerthresh)) {
            # If the power never gets high enough - take the max and alert with a message
            powerthresh <- sft$fitIndices[sft$fitIndices[,2] == max(sft$fitIndices[,2]),1]
        }
    } else  {
        powerthresh = powerparam
    }
    
    ## Plot the powerplot
    powerplot <- scatter_plotter(sft$fitIndices[,c(1,2)], colorvar=data.frame(PowerEstimate = sft$fitIndices[,c(1)] == powerthresh), 
                                 shapevar = NULL, sizevar = NULL, datalabels = NULL,
                                 labsparam = list(title = "power estimation plot", x = "Power", y = "Fit", color = "PowerEstimate"), plotstats = FALSE)
    # pdf(paste0(outfilepath, "powerplot.pdf"))
    # print(powerplot)
    # junk <- dev.off()
    
    
    ## Creating modules is a single call
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    cor <- WGCNA::cor
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    net = blockwiseModules(counttable, power = powerthresh, 
                           maxBlockSize = (ncol(counttable)+1),
                           # blocks = 1,
                           TOMType = "unsigned", minModuleSize = minModuleSizeParam,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = FALSE,
                           verbose = 3)
    

    
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    cor<-stats::cor
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    # The number and size of the modules: NOTE - 0 is for genes outside ALL modules
    moduleLabels <- net$colors
    moduleColors = labels2colors(moduleLabels)
    genestocolors <- cbind.data.frame(moduleLabels, moduleColors, stringsAsFactors = FALSE)
    table(genestocolors[,1])
    table(genestocolors[,2])

    if (!is.null(savedendroplotfile)) {
        consTree <- net$dendrograms
        pdf(savedendroplotfile)
        plotDendroAndColors(consTree[[1]], moduleColors,
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Consensus gene dendrogram and module colors")
        junk <- dev.off()
    }

    
    
    
    
    
    
    return(list(powerplot = powerplot, genestocolors = genestocolors))
}


# annotate_WGCNA_modules <- function(counttable, genestocolors, genesetparam = "C5") {
annotate_WGCNA_modules <- function(genestocolors, genesetparam = "C5", padjcutoffval = 0.1) {
    # moduleColors <- genestocolors[,2]
    
    # Create eigengenes
    # eigengenes = orderMEs(moduleEigengenes(counttable, moduleColors)$eigengenes)
    
    ## Align modules with known GO terms
    # genestocolors
    # modNames = substring(names(eigengenes), 3)
    modNames = unique(genestocolors[,2])
    outenrichmentlist <- list()
    for (modulecolornum in seq_len(length(modNames))) {
        
        ## Pick the module
        module = modNames[modulecolornum]
        GOI <- data.frame(GOI = rownames(genestocolors[genestocolors[,2] == module,2,drop=FALSE]), 
                          row.names = rownames(genestocolors[genestocolors[,2] == module,2,drop=FALSE]))
        
        out1 <- hypergeo_genetest(DEseqtable = GOI, statcutoffparam = NULL, genesetparam = c(genesetparam), speciesparam = "Homo sapiens")[[1]]
        if (nrow(out1) == 0) {
            out1[1,] <- NA
            numpathways <- 1
        } else {
            # Only keep a certain number of hits - either a min of 10, max of 100, p.adjust < 0.1
            padjcutoffval <- padjcutoffval
            # numpathways <- min(100, max(10, nrow(out1[out1[,"p.adjust"] < padjcutoffval])))
            numpathways <- max(10, nrow(out1[out1[,"p.adjust"] < padjcutoffval]))
        }

        # Create the outtable
        outtab <- cbind(module = module, size = nrow(GOI), out1[,c("pvalue", "p.adjust", "GeneRatio")], ont = genesetparam, out1[,"ID", drop = FALSE])
        outtab <- outtab[1:numpathways,]
        
        # Save it an outlist to bind together at the end
        outenrichmentlist[[modulecolornum]] <- outtab
        names(outenrichmentlist)[modulecolornum] <- module
        
    }
    moduleenrichmenttab <- do.call(rbind, outenrichmentlist)
    # write.table(moduleenrichmenttab, paste0(outfilepath, "module_enrichment_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    return(moduleenrichmenttab)
}


# WGNCA_module_feature_correlation <- function(counttable, moduleColors, metatable, outfilepath) {
WGNCA_module_feature_correlation <- function(counttable, moduleColors, metatable) {
    ## Create summary heatmap of the correlation results
    ## PART 3 - RELATING MODULES TO EXISTING INFORMATION AND IDENTIFYING KEY GENES
    # the eigengene is the summary profile for each module, we can then correlate this info with measured clinical traits to look for associations
    print(paste0("Check for same order of metatable and counttable: ", identical(rownames(counttable), rownames(metatable))))
    
    # Create eigengenes
    eigengenes = orderMEs(moduleEigengenes(counttable, moduleColors)$eigengenes)
    
    if (sum(unlist(lapply(metatable, function(x) is.character(x)|is.factor(x)))) > 0) {
        out1 <- apply(make_comparison_columns(metatable[,unlist(lapply(metatable, function(x) is.character(x)|is.factor(x))),drop=FALSE]),
                      2, function(x) as.numeric(factor(x)))
        rownames(out1) <- rownames(metatable)
        out2 <- metatable[,!colnames(metatable) %in% colnames(metatable[,unlist(lapply(metatable, function(x) is.character(x)|is.factor(x))),drop=FALSE])]
        WGCNA_meta <- cbind(out1, out2)
    } else {
        WGCNA_meta <- metatable
    }
    
    ## Calculate the module feature correlation table
    moduleTraitCor = cor(eigengenes, WGCNA_meta, use = "p", method = "spearman")
    # write.table(moduleTraitCor, paste0(outfilepath, "module_trait_R_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(counttable))
    # write.table(moduleTraitPvalue, paste0(outfilepath, "module_trait_p_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    
    ## Heatmap
    # DOESNT MATTER, JUST A HEATMAP OFF THE CORRELATION VALUES - WHILE USEFUL, I CAN DO SOMETHING ON MY OWN
    
    # pdf(paste0(outfilepath, "correlation_heatmap.pdf"), 12, 12)
    # # sizeGrWindow(20,12)
    # # Will display correlations and their p-values
    # textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
    #                     signif(moduleTraitPvalue, 1), ")", sep = "")
    # dim(textMatrix) = dim(moduleTraitCor)
    # # par(mar = c(6, 8.5, 3, 3))
    # par(mar = c(13, 10, 3, 3))
    # # Display the correlation values within a heatmap plot
    # labeledHeatmap(Matrix = moduleTraitCor,
    #                xLabels = names(WGCNA_meta),
    #                yLabels = names(eigengenes),
    #                ySymbols = names(eigengenes),
    #                colorLabels = FALSE,
    #                # colors = blueWhiteRed(50),
    #                colors = NULL,
    #                textMatrix = textMatrix,
    #                setStdMargins = FALSE,
    #                cex.text = 0.5,
    #                zlim = c(-1,1),
    #                main = paste("Module-trait relationships"))
    # junk <- dev.off()
    
    return(list(moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue))
}














run_WGCNA <- function(counttable, metatable, powerparam = NULL, outfilepath){

    ## Create outfolder in case it wasnt
    dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)

    powers = c(1:20)
    sft = pickSoftThreshold(counttable, powerVector = powers, verbose = 5)
    ## Determine the power needed
    if (is.null(powerparam)) {
        powerthresh = sft$powerEstimate
        if (is.na(powerthresh)) {
            # If the power never gets high enough - take the max and alert with a message
            powerthresh <- sft$fitIndices[sft$fitIndices[,2] == max(sft$fitIndices[,2]),1]
        }
    } else  {
        powerthresh = powerparam
    }
    
    ## Plot the powerplot
    powerplot <- scatter_plotter(sft$fitIndices[,c(1,2)], colorvar=data.frame(PowerEstimate = sft$fitIndices[,c(1)] == powerthresh), 
                    shapevar = NULL, sizevar = NULL, datalabels = NULL,
                    labsparam = list(title = "power estimation plot", x = "Power", y = "Fit", color = "PowerEstimate"), plotstats = FALSE)
    pdf(paste0(outfilepath, "powerplot.pdf"))
    print(powerplot)
    junk <- dev.off()
    
    
    ## Creating modules is a single call
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    cor <- WGCNA::cor
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    net = blockwiseModules(counttable, power = powerthresh,
                           TOMType = "unsigned", minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = FALSE,
                           verbose = 3)
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    cor<-stats::cor
    ### THIS IS EFFING AMAZING - but theres a bug in WGNCA - you have to assign and then revert the cor function
    # The number and size of the modules: NOTE - 0 is for genes outside ALL modules
    moduleLabels <- net$colors
    moduleColors = labels2colors(moduleLabels)
    genestocolors <- cbind.data.frame(moduleLabels, moduleColors, stringsAsFactors = FALSE)
    table(genestocolors[,1])
    table(genestocolors[,2])
    # geneTree = net$dendrograms[[1]]
    eigengenes = orderMEs(moduleEigengenes(counttable, moduleColors)$eigengenes)
    
    write.table(genestocolors, paste0(outfilepath, "genestocolors.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    
    ## Align modules with known GO terms
    # genestocolors
    modNames = substring(names(eigengenes), 3)
    outenrichmentlist <- list()
    for (modulecolornum in seq_len(length(modNames))) {
        
        ## Pick the module
        module = modNames[modulecolornum]
        GOI <- data.frame(GOI = rownames(genestocolors[genestocolors[,2] == module,2,drop=FALSE]), 
                          row.names = rownames(genestocolors[genestocolors[,2] == module,2,drop=FALSE]))
        
        out1 <- hypergeo_genetest(GOI, statcutoffparam = NULL, genesetparam = c("C5"), speciesparam = "Homo sapiens")[[1]]
        
        # Only keep a certain number of hits - either a min of 10, max of 100, p.adjust < 0.1
        padjcutoffval <- 0.1
        # numpathways <- min(100, max(10, nrow(out1[out1[,"p.adjust"] < padjcutoffval])))
        numpathways <- max(10, nrow(out1[out1[,"p.adjust"] < padjcutoffval]))
        
        # Create the outtable
        outtab <- cbind(module = module, size = nrow(GOI), out1[,c("pvalue", "p.adjust", "GeneRatio")], ont = "C5", out1[,"ID", drop = FALSE])
        outtab <- outtab[1:numpathways,]
        
        # Save it an outlist to bind together at the end
        outenrichmentlist[[modulecolornum]] <- outtab
        names(outenrichmentlist)[modulecolornum] <- module
        
    }
    moduleenrichmenttab <- do.call(rbind, outenrichmentlist)
    write.table(moduleenrichmenttab, paste0(outfilepath, "module_enrichment_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    
    
    moduleenrichmenttab[grepl("PLATELET", moduleenrichmenttab[,"ID"]),]
    
    
    
    
    
    ## Create summary heatmap of the correlation results
    ## PART 3 - RELATING MODULES TO EXISTING INFORMATION AND IDENTIFYING KEY GENES
    # the eigengene is the summary profile for each module, we can then correlate this info with measured clinical traits to look for associations
    # Define numbers of genes and samples
    # nGenes = ncol(counttable)
    # nSamples = nrow(counttable)
    # Recalculate MEs with color labels
    # MEs0 = moduleEigengenes(counttable, moduleColors)$eigengenes
    # MEs = orderMEs(MEs0)
    
    if (sum(unlist(lapply(metatable, function(x) is.character(x)|is.factor(x)))) > 0) {
        out1 <- apply(make_comparison_columns(metatable[,unlist(lapply(metatable, function(x) is.character(x)|is.factor(x)))]),
                      2, function(x) as.numeric(factor(x)))
        rownames(out1) <- rownames(metatable)
        out2 <- metatable[,!colnames(metatable) %in% colnames(metatable[,unlist(lapply(metatable, function(x) is.character(x)|is.factor(x)))])]
        WGCNA_meta <- cbind(out1, out2)
    } else {
        WGCNA_meta <- metatable
    }
    
    
    
    ## Calculate the module feature correlation table
    moduleTraitCor = cor(eigengenes, WGCNA_meta, use = "p", method = "spearman")
    write.table(moduleTraitCor, paste0(outfilepath, "module_trait_R_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(counttable))
    write.table(moduleTraitPvalue, paste0(outfilepath, "module_trait_p_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    
    ## Heatmap
    # DOESNT MATTER, JUST A HEATMAP OFF THE CORRELATION VALUES - WHILE USEFUL, I CAN DO SOMETHING ON MY OWN
    
    pdf(paste0(outfilepath, "correlation_heatmap.pdf"), 12, 12)
    # sizeGrWindow(20,12)
    # Will display correlations and their p-values
    textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    # par(mar = c(6, 8.5, 3, 3))
    par(mar = c(13, 10, 3, 3))
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = names(WGCNA_meta),
                   yLabels = names(eigengenes),
                   ySymbols = names(eigengenes),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    junk <- dev.off()


    
}


## CUSTOM WGCNA HEATMAP
WGCNA_custom_heatmap <- function(moduleTraitPvalue, moduleTraitCor) {
    maptab = -log10(as.matrix(moduleTraitPvalue)) * ifelse(as.matrix(moduleTraitCor) > 0, 1, -1)
    # rownames(maptab) <- str_wrap(gsub("_", " ", gsub("^GO|^KEGG|^HALLMARK", "", rownames(maptab))), 20)
    colannotationlist = NULL
    colmetatable = NULL
    
    ## Build Annotations from our metatable
    hatop = NULL
    if (!is.null(colannotationlist) & !is.null(colmetatable)) {
        temp1 <- vector("list", length(colannotationlist))
        names(temp1) = names(colannotationlist)
        annotlegendlist = lapply(temp1, function(x) x[[1]] = 
                                     list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
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
                                  show_legend = TRUE,
                                  annotation_legend_param = annotlegendlist)
    }
    # Define the colormap
    # heatmapcolorparam = colorRamp2(
    #     # c(floor(min(maptab)), -0.99999, 0, 0.99999, ceiling(max(maptab))), c("blue", "white", "white", "white", "red"))
    #     c(floor(min(maptab)), -log10(0.2), 0, log10(0.2), ceiling(max(maptab))), c("blue", "white", "white", "white", "red"))
    
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
    heatmapcolorparam = colorRamp2(
        c(lowvalue, midvalue, highvalue), c(lowcolor, midcolor, highcolor))
    
    ht1 = Heatmap(maptab, 
                  border = TRUE,
                  rect_gp = gpar(col = "black", lwd = 1),
                  col = heatmapcolorparam,    ## Define the color scale for the heatmap
                  row_title = "Factor",                                       ## Name the rows
                  column_title = "Modules",                                  ## Name the columns
                  
                  cluster_columns = FALSE,                         ## Cluster the columns or leave as is
                  cluster_rows = FALSE,                            ## Cluster the rows or leave as is
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
                  
                  show_column_names = TRUE  ,                               ## Show the Column Names
                  column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
                  show_row_names = TRUE,                                   ## Show the row names
                  row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
                  row_names_gp = gpar(fontsize=6),
                  
                  show_row_dend = TRUE, #nrow(maptab) <=500,                                     ## Show the dendrogram on the rows
                  show_column_dend = TRUE,                                   ## Show the dendrogram on the columns
                  
                  heatmap_legend_param = list(
                      title = "-log10(pvalue)",
                      legend_height = unit(2.5, "cm"),
                      title_gp = gpar(fontsize = 8, fontface = "bold")),
                  # width = unit(4,"cm"),
                  # height = unit(20,"cm"),
                  # height = unit(min((nrow(maptab)/2), 20),"cm"),
                  # width = unit(min(ncol(maptab), 30),"cm"),
                  # height = unit(min((nrow(maptab)*1.8), 30),"cm"),
                  height = unit(max((nrow(maptab)*1.5), 20),"cm"),
                  width = unit(min(ncol(maptab), 30),"cm"),
                  
                  # top_annotation = hatop,
                  
                  ## Adds in the cor value from the maptab
                  # cell_fun = function(j, i, x, y, width, height, fill) {
                  #     grid.text(paste0(sprintf("%.2f", t(moduleTraitCor)[,colnames(maptab),drop=FALSE][i, j]), "\n", "(",
                  #                      sprintf("%.2f", t(moduleTraitPvalue)[,colnames(maptab),drop=FALSE][i, j]), ")"),
                  #               x, y, gp = gpar(fontsize = 8))
                  # }
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(paste0(sprintf("%.2f", moduleTraitCor[i, j]), "\n", "(",
                                       sprintf("%.2f", moduleTraitPvalue[i, j]), ")"),
                                x, y, gp = gpar(fontsize = 8))
                  }
                  
    )
    
    return(ht1)
    # pdf(paste0(WGCNA_outfilepath, "feature_module_corr_hm_custom.pdf"), 6, 14)
    # draw(ht1)
    # junk <- dev.off()
}





## Need to have selected pathways and color order preset, will not do in function
WGCNA_annotation_dotplot <- function(WGCNA_annotation_table, padjcutoffparam = 0.05, pathwayselect = NULL, cleannames = FALSE) {
    
    ## I need to think about how to add a side annotation or the genesets that are a part of each module
    modulesel <- na.omit(WGCNA_annotation_table[, c("module", "p.adjust", "ID")])
    if (!is.null(padjcutoffparam)) {
        modulesel <- modulesel[modulesel[,"p.adjust"] < padjcutoffparam,]
        padj_plot_upperlimit <- padjcutoffparam
    } else {
        padj_plot_upperlimit <- 1
    }
    if (!is.null(pathwayselect)) {modulesel <- modulesel[modulesel[,"ID"] %in% pathwayselect,]}
    
    # add on the other colors with blanks
    if (sum(!unique(WGCNA_annotation_table[,"module"]) %in% modulesel[,"module"]) != 0) {
        blankmoduletab <- data.frame(module = unique(WGCNA_annotation_table[,"module"])[!unique(WGCNA_annotation_table[,"module"]) %in% modulesel[,"module"]],
                                     p.adjust = NA, ID = modulesel[1,"ID"])
    } else {blankmoduletab <- NULL}
    dotplottab <- rbind(modulesel, blankmoduletab)
    
    ## Reorder to match the order of eigengenes
    dotplottab <- dotplottab[order(dotplottab[,"module"], dotplottab[,"p.adjust"]),]
    ## Clean names if turned on
    if (cleannames == TRUE){
        dotplottab[,3] <- gsub("HALLMARK_|GO_|HP_", "", dotplottab[,3])
        dotplottab[,3] <- unlist(lapply(gsub("_", " ", gsub("HALLMARK_", "", dotplottab[,3])), function(x) simpleCap(tolower(x))))
        dotplottab[,3] <- factor(dotplottab[,3], levels = unique(dotplottab[,3]))
    } else {
        dotplottab[,3] <- factor(dotplottab[,3], levels = unique(dotplottab[,3]))
    }
    
    labsparam <- list(title = "module dot plot", x = "modules", y = "GO terms", color = "module", size = "adjusted p value")
    colorvec <- as.character(unique(dotplottab[,1]))
    names(colorvec) <- colorvec
    
    pout <- ggplot(dotplottab, aes(x=dotplottab[,1], y=dotplottab[,3], size=dotplottab[,2], color=as.character(dotplottab[,1])))
    # pout <- ggplot(modulesel, aes(x=modulesel[,1], y=modulesel[,3], size=modulesel[,2], color=modulesel[,2]))
    pout <- pout + geom_point(alpha = 0.8, shape = 16)
    pout <- pout + scale_size_continuous(name = "adjusted p value", range = c(8, 2), limits = (c(0,padj_plot_upperlimit)))
    pout <- pout + scale_color_manual(values = colorvec, guide = FALSE)
    pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, color = labsparam$color, shape = labsparam$shape, size = labsparam$size)
    pout <- pout + theme_bw()
    pout <- pout + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))
    
    return(list(dotplot = pout, dotplottab = dotplottab))

}




