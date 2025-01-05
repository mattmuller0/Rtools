###########################################################################
#
#                            proteomics_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: proteomics_functions

#======================== LIBRARIES ========================#
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
suppressMessages(library(OlinkAnalyze))
suppressMessages(library(clusterProfiler))
suppressMessages(library(AnnotationDbi))

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')

message("
     | \          .-'               |    \     .- |__
     |  -,-.-.   '-. _  \ .\ .  .   |     || | `. |
     |  |  | |     || `\   |  \ |   |  `._/`-'_.' |
     '          '-' `._/     .  |   |             `-'
                       |      `-'                 _  .
      _  _  .  \ |__   |     _     .-.  _     _ .' | |
     /  | `||,-. |     |--. /-'   |  __| `\ /' ||  | |
     ._ `-'`|  | |     |   |\ ,   |   |`._/ \ .' \ | '
               ` '-    |`-'__      `-'               o
                      ,-'`` .'
       /`.           /   .'`    .-;   .-.
      /   \          |   \   _,' /-,_  \ `.
     / _ .-'    _.---:.    `'    _.-'_.'   )
     '` \\   ,;'               .  \ `     ,'
         ')    `7_,             `'-;   .-'-.._
        //   ,;-'`        .-''. ,-.        _.f`
      .'/   .-'`\ _      .-'`\ \ _ `'-       >
     / /   /    /` \   ,'     ; / `'.      ,-''-.
    |  |  |     \   \ /       |;     \   ,'      `.
    \   '.;      `.  \ -.     ||       _/         |
_mx__`.___\\    / L\  \  \ _O_;|O__ .'`          /________
            `.  |  `'-.(  | _..--''/         _.-`
              `-.`.__.p \ -'       \ /  ,  ;`          /.
           '.----'L__p\ '.          v|  |   \ `'-._   /  \
             `.   '.__b\  \       _~  v-'`v-' ,\'  `.'   |
               `.       `<~|   .-'  ,' _ _.'`   _.'     /
                 `.       `' '_,-''  ,.-''`-  _.'      .'
                   `.       \   _,,='      ,-'     _.-'
                     `.         _  __ _.-'     _.-'
                       `.      `. `  `     _.-'
                         `.      \     _.-'
                           `.     |_.-'
                             `-._.'
")

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
        { if (byPanel) dplyr::select(., SampleID, Outlier, Panel) else dplyr::select(., SampleID, Outlier) }
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
        { if (byPanel) dplyr::select(., SampleID, Outlier, Panel) else dplyr::select(., SampleID, Outlier) }
    write.csv(outliers, glue('{outdir}/umap_outliers.csv'), row.names = F)

    out <- list(umap = umap, outliers = outliers)
    return(out)
}

# Function test level of detection (LOD) for each protein
olink_lod_qc <- function(data, outdir, plot = TRUE) {
    dir.create(outdir, showWarnings = F)
    # get the lod
    lod <- data %>% dplyr::mutate(LOD_QC = NPX > LOD)
    write.csv(lod, glue('{outdir}/lod_qc.csv'), row.names = F)
    
    # plot the NPX facet wrappped by Assay and colored by LOD_QC
    if (plot) {
        p <- ggplot(lod, aes(x = NPX, fill = LOD_QC)) +
            geom_histogram(bins = 100) +
            facet_wrap(~Assay, scales = 'free') +
            theme_bw() +
            theme(legend.position = 'bottom')
        ggsave(glue('{outdir}/npx_hist.pdf'), p, width = 45, height = 25)
    }

    # proteins that fail the LOD in more than 75% of samples
    # outliers <- lod %>%
    #     dplyr::group_by(Assay) %>%
    #     dplyr::summarise(n_pass = sum(LOD_QC) / n()) %>%
    #     dplyr::filter(n_pass < 0.75) 
    outliers <- lod %>% dplyr::filter(MissingFreq > 0.25)

    # return the proteins that fail the LOD
    out <- list(lod = lod, outliers = outliers)
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

# Function to do olink filtering
# Arguments:
#   data: data frame of olink data
olink_filtering <- function(data, outdir, pca_args = list(), umap_args = list(), lod_args = list()) {
    dir.create(outdir, showWarnings = F)
    # get the sample info
    message('Getting sample info')
    sample_info <- do.call(olink_info, list(data))
    writeLines(glue("{names(sample_info)}: {sample_info}\n"), file.path(outdir, "sampleIDs.txt"))

    # get the pca outliers
    message('Running PCA')
    pca_outliers <- do.call(olink_pca_outliers, c(list(data = data, outdir = file.path(outdir, 'pca')), pca_args))
    saveRDS(pca_outliers, glue('{outdir}/pca/pca_outliers.rds'))

    # get the umap outliers
    message('Running UMAP')
    umap_outliers <- do.call(olink_umap_outliers, c(list(data = data, outdir = file.path(outdir, 'umap')), umap_args))
    saveRDS(umap_outliers, glue('{outdir}/umap/umap_outliers.rds'))

    # get the lod qc
    message('Running LOD')
    lod_qc <- do.call(olink_lod_qc, c(list(data = data, outdir = file.path(outdir, 'lod')), lod_args))
    saveRDS(lod_qc, glue('{outdir}/lod/lod_qc.rds'))

    # get the outliers
    message('Getting outliers')
    qc_outliers <- data %>% filter(QC_Warning == "WARN") %>% pull(SampleID) %>% unique()
    pca_outliers <- pca_outliers$outliers %>% pull(SampleID)
    umap_outliers <- umap_outliers$outliers %>% filter(Outlier == 1) %>% pull(SampleID)
    outliers <- unique(c(qc_outliers, pca_outliers, umap_outliers))
    write.csv(outliers, glue('{outdir}/sample_outliers.csv'), row.names = F)

    # proteins that are below the LOD
    message('Getting proteins below the LOD')
    lod_outliers <- lod_qc$outliers %>% pull(Assay) %>% unique()
    write.csv(lod_outliers, glue('{outdir}/protein_outliers.csv'), row.names = F)

    # remove the outliers
    data %<>% filter(!grepl(paste(outliers, collapse = "|"), SampleID))
    # remove the proteins that are below the LOD
    data %<>% filter(!grepl(paste(lod_outliers, collapse = "|"), Assay))
    write.csv(data, glue('{outdir}/npx_data.csv'), row.names = F)

    # get the count table
    count_table <- olink_count_table(data)
    write.csv(count_table, glue('{outdir}/count_table.csv'))

    # return the data
    message('Done filtering')
    out <- list(data = data, count_table = count_table, outliers = outliers, lod_outliers = lod_outliers)

    return(out)
}

# Function to run differential expression and pathway analysis
# Arguments:
#   data: data frame of olink data
#   condition: condition column
#   outdir: output directory
#   de_fxn: differential expression function
#   gsea_fxn: pathway analysis function
#   de_args: arguments to pass to differential expression function
#   gsea_args: arguments to pass to pathway analysis function
# Outputs:
#   list of differential expression and pathway analysis results
olink_analysis <- function(
    data, 
    conditions, 
    outdir, 
    de_fxn = olink_wilcox, 
    gsea_fxn = olink_pathway_enrichment,
    de_args = list(),
    volcano_args = list(),
    gsea_args = list()

    ) {
    dir.create(outdir, showWarnings = F)
    
    # make sure the conditions are in the data
    if (!all(conditions %in% colnames(data))) {
        stop('Conditions not in data')
    }

    de_res <- list()
    gsea_res <- list()
    res <- list()
    message("Looping over: ", paste0(conditions, collapse = ", "))
    for (condition in conditions) {
        dir.create(glue('{outdir}/{condition}'), showWarnings = F)
        # get the differential expression
        message('Running differential expression')
        de <- do.call(de_fxn, c(list(data, condition), de_args))
        p <- do.call(olink_volcano_plot, list(de))
        ggsave(glue('{outdir}/{condition}/volcano.pdf'), p, width = 6, height = 6)
        write.csv(de, glue('{outdir}/{condition}/de_results.csv'), row.names = F)
        de_res[[condition]] <- de
        r <- summarize_experiment(de, logFC_column = "estimate", pvalue_column = "p.value", padj_column = "Adjusted_pval")
        r$condition <- condition
        res[[condition]] <- r

        # get the pathway analysis
        message('Running pathway analysis')
        gsea <- do.call(olink_pathway_enrichment, c(list(data, de), gsea_args))
        write.csv(gsea, glue('{outdir}/{condition}/gsea_results.csv'), row.names = F)
        p <- do.call(olink_pathway_visualization, list(gsea))
        ggsave(glue('{outdir}/{condition}/gsea.pdf'), p, width = 12, height = 8)
        p <- do.call(olink_pathway_heatmap, list(gsea, de))
        ggsave(glue('{outdir}/{condition}/gsea_heatmap.pdf'), p, width = 12, height = 8)
        gsea_res[[condition]] <- gsea
    }

    res <- do.call(rbind, res)
    write.csv(res, glue('{outdir}/experiment_summary.csv'), row.names = F)

    # return the results
    message('Done with analysis')
    out <- list(de = de_res, gsea = gsea_res)
    return(out)
}

# Function to calculate odds ratios for olink data
# Arguments:
#   data: data frame of olink data
#   proteins: vector of proteins to calculate odds ratios for
#   events: vector of events to calculate odds ratios for
#   adjustments: vector of variables to adjust for
# Outputs:
#   list of odds ratios for each event
olink_odds_ratios <- function(data, proteins, events, adjustments = NULL) {
    # map over the events and proteins to get the odds ratios
    or.df <- purrr::map(
        .x = events,
        .f = function(event) {
            out <- purrr::map(
                .x = proteins,
                .f = function(protein) {
                    # get the counts for the event and protein
                    dat <- data %>%
                        dplyr::select(protein, event, any_of(adjustments)) %>%
                        drop_na()
                    # make a formula to adjust for
                    if (is.null(adjustments)) {
                        form <- paste0(event, ' ~ ', protein)
                    } else {
                        form <- paste0(event, ' ~ ', protein, ' + ', paste(adjustments, collapse = ' + '))
                    }
                    # get the odds ratio
                    odds_ratio <- glm(
                        formula = form,
                        data = dat,
                        family = binomial(link = 'logit')
                    ) %>%
                        broom::tidy() %>%
                        filter(term == protein) %>%
                        mutate(
                            odds.ratio = exp(estimate),
                            or.ci.lower = exp(estimate - 1.96 * std.error),
                            or.ci.upper = exp(estimate + 1.96 * std.error),
                            formula = form
                        ) %>%
                        dplyr::select(term, odds.ratio, or.ci.lower, or.ci.upper, estimate, std.error, p.value, formula)
                }
            )
            out <- do.call(rbind, out)
            out$p.adj <- p.adjust(out$p.value, method = 'BH')
            out
        }
    )
    names(or.df) <- events
    return(or.df)
}
