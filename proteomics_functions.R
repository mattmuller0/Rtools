###########################################################################
#
#                            proteomics_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: proteomics_functions

#======================== LIBRARIES ========================#
library(lintr) #nolint
library(httpgd) #nolint
library(languageserver) #nolint
library(devtools)
library(tidyverse)
library(glue)
library(OlinkAnalyze)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')

#======================== CODE ========================#

# Function to get olink sample info
olink_info <- function(data) {
    # get the sample IDs
    sampleIDs <- unique(data$SampleID)

    # get the plates
    plates <- unique(data$PlateID)

    # get the proteins
    proteins <- unique(data$Assay)

    # get the protein groups
    panels <- unique(data$Panel)

    # write out the sample IDs, plates, proteins, and protein groups to a file
    out <- list(
        sampleIDs = sampleIDs,
        plates = plates,
        proteins = proteins,
        panels = panels
        )
    return(out)
}

# Function to get olink PCA outliers and plot
# Arguments:
#   data: data frame of olink data
#   outdir: output directory
#   outlierDefX: outlier definition for x axis
#   outlierDefY: outlier definition for y axis
#   byPanel: boolean to plot by panel
#   ...: additional arguments to pass to olink_pca_plot
# Outputs:
#   plot of pca outliers
olink_pca_outliers <- function(data, outdir, outlierDefX = 2.5, outlierDefY = 4, byPanel = TRUE, ...) {
    dir.create(outdir, showWarnings = F)
    # get the pca
    pca <- olink_pca_plot(data, outlierDefX = outlierDefX, outlierDefY = outlierDefY, byPanel = byPanel, ...)
    lapply(names(pca), function(x) {ggsave(glue('{outdir}/{x}_pca.pdf'), pca[[x]])})

    # outliers
    outliers <- lapply(pca, function(x) {x$data}) %>%
        bind_rows() %>%
        filter(Outlier == 1) %>% 
        dplyr::select(SampleID, Outlier, Panel)
    write.csv(outliers, glue('{outdir}/pca_outliers.csv'), row.names = F)

    out <- list(pca = pca, outliers = outliers)
    return(out)
}

# Function to get olink UMAP outliers and plot
# Arguments:
#   data: data frame of olink data
#   outdir: output directory
#   outlierDefX: outlier definition for x axis
#   outlierDefY: outlier definition for y axis
#   byPanel: boolean to plot by panel
#   ...: additional arguments to pass to olink_pca_plot
# Outputs:
#   plot of pca outliers
olink_umap_outliers <- function(data, outdir, outlierDefX = 2.5, outlierDefY = 4, byPanel = TRUE, ...) {
    dir.create(outdir, showWarnings = F)
    # get the pca
    umap <- olink_umap_plot(data, outlierDefX = outlierDefX, outlierDefY = outlierDefY, byPanel = byPanel, ...)
    lapply(names(umap), function(x) {ggsave(glue('{outdir}/{x}_umap.pdf'), umap[[x]])})

    # outliers
    outliers <- lapply(umap, function(x) {x$data}) %>%
        bind_rows() %>%
        filter(Outlier == 1) %>% 
        dplyr::select(SampleID, Outlier, Panel)
    write.csv(outliers, glue('{outdir}/umap_outliers.csv'), row.names = F)

    out <- list(umap = umap, outliers = outliers)
    return(out)
}

# Function test level of detection (LOD) for each protein
olink_lod_qc <- function(data, outdir) {
    dir.create(outdir, showWarnings = F)
    # get the lod
    lod <- data %>% dplyr::mutate(LOD_QC = NPX > LOD)
    
    # plot the NPX facet wrappped by Assay and colored by LOD_QC
    p <- ggplot(lod, aes(x = NPX, fill = LOD_QC)) +
        geom_histogram(bins = 100) +
        facet_wrap(~Assay, scales = 'free') +
        theme_bw() +
        theme(legend.position = 'bottom')
    ggsave(glue('{outdir}/npx_hist.pdf'), p, width = 45, height = 25)

    # get the number of proteins that pass the LOD
    n_pass <- lod %>% dplyr::filter(LOD_QC == TRUE) %>% dplyr::select(Assay) %>% unique() %>% nrow()
    message(glue('Number of proteins that pass the LOD: {n_pass}'))
    # get the number of proteins that fail the LOD
    n_fail <- lod %>% dplyr::filter(LOD_QC == FALSE) %>% dplyr::select(Assay) %>% unique() %>% nrow()
    message(glue('Number of proteins that fail the LOD: {n_fail}'))

    # write out the lod data
    write.csv(lod, glue('{outdir}/lod_qc.csv'), row.names = F)

    # return the proteins that fail the LOD
    out <- list(lod = lod, n_pass = n_pass, n_fail = n_fail, qc_plot = p)
    return(out)
}


# Function to convert olink data to count table
# Arguments:
#   data: data frame of olink data
#   sampleID: sample ID column
#   assay: assay column
#   value: value column
# Outputs:
#   count table
olink_count_table <- function(data, sampleID = 'SampleID', assay = 'Assay', value = 'NPX') {
    # get the count table
    count_table <- data %>%
        # get the sample ID, assay, and value
        dplyr::select(sampleID, assay, value) %>%
        # pivot to wide
        tidyr::pivot_wider(names_from = assay, values_from = value) %>%
        # remove the sample ID
        column_to_rownames(sampleID)
    return(count_table)
}