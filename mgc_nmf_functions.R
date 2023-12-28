################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## NMF Functions

# https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
# https://cran.r-project.org/web/packages/NMF/vignettes/heatmaps.pdf

#### NOTE - DUE TO THE BUG OF NMF PACKAGES, THIS MUST BE RUN BYITSELF DUE TO ITS USE OF SEEDING AND CONFLICT WITH OTHER PACKAGES
#### SO YOU MUST RUN THIS ON ITS OWN, NOT LOADED WITH A BUNCH OF OTHER THINGS

library("NMF")

## Step 1 - clean data as needed
# nmf_format_data <- function(){}

## Step 1 - clean data as needed
nmf_ranktest <- function(incounttable, ranksparam = 2:8, nrunparam = 50, seedparam = 123456, plotconsensusmapparam = NULL){
    
    ##  To determine the optimal factorization rank we calculated two metrics for each k:
    # 1) cophenetic correlation coefficient measuring how well the intrinsic structure of the data was recapitulated after clustering and
    # 2) the dispersion coefficient of the consensus matrix as defined in Kim and Park (2007) measuring the reproducibility of the clustering across 50 iterations. 
    # 3) The optimal k was defined as the maximum of the product of both metrics for cluster numbers between k = 3 and 8 (FiguresS3C and S4B).
    
    ranktest_out <- nmf(incounttable, ranksparam, nrun=nrunparam, seed=seedparam, .opt='v4')
    # saveRDS(object = ranktest_out, file = paste0(outfilepathnmf, "ranktest_object.rds"))
    rankmetricplot <- plot(ranktest_out) ### SAVE THIS OUT
    # pdf(paste0(outfilepathnmf, "rankmetric_plot.pdf"), useDingbats = FALSE)
    # print(rankmetricplot)
    # dev.off()

    ## So they multiplied the cophenetic and dispersion value together - and then took the max for k 3-8 (not 2 because 2 will always be highest!)
    ranktest_measures <- ranktest_out$measures
    ranktest_measures[,"ranktest_metric"] <- ranktest_measures[,"cophenetic"] * ranktest_measures[,"dispersion"]
    suggested_rank <- ranktest_measures[ranktest_measures[,"ranktest_metric"] == max(ranktest_measures[ranktest_measures[,"rank"] > 2,"ranktest_metric"]), "rank"] ## SAVE THIS OUT
    write.table(ranktest_measures, paste0(outfilepathnmf, "rankmetric_table.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
    
    ## estim.r is a list with measures, consensus, fit
    ranktest_out$measures # data.frame with quality measures in column, ranks per row
    ranktest_out$consensus # list of consensus matrices for each rank
    ranktest_out$fit # the output metric for each nmf run
    
    ## Save this out too
    if (!is.null(plotconsensusmapparam)) {
        dir.create(plotconsensusmapparam, showWarnings = FALSE, recursive = TRUE)
        pdf(paste0(plotconsensusmapparam, "consensus_matrix_plot.pdf"), width = 20, height = 12, useDingbats = FALSE)
        consensusmap(ranktest_out, labCol=NA, labRow=NA)
        dev.off()
    }
    
    ## Outlist
    return(list(ranktestobject = ranktest_out,
                rankmetricplot = rankmetricplot,
                rankmetrictable = ranktest_measures,
                suggestedrank = suggested_rank))
    
}


run_nmf <- function(incounttable, ranksparam, methodparam = NULL, nrunparam = 200) {
    
    ## Run method
    inmethod = ifelse(is.null(methodparam), "brunet", methodparam)
    nmfout_object <- nmf(incounttable, rank = ranksparam, method = inmethod, nrun=nrunparam, .opt='v4')
    # saveRDS(object = nmfout, file = paste0(outfilepathnmf, "nmf_object.rds"))
    nmf_fitted_matrix <- fitted(nmfout_object)
    summary(nmfout_object)
    
    #
    nmf_w <- basis(nmfout_object) # feature x rank - weights for each feature (p x k, p is feature number)()
    nmf_h <- coef(nmfout_object) # rank x samples - rank weights for each sample (k x n, n is feature number)()
    ## From nmf_h - H is used to assign samples to clusters by choosing the k with max score in each col (so max k value for each sample)
    ## Each sample gets a cluster membership score - by calculating the maximal fractional score of the corresponding column.
    nmf_h_membership_scores <- t(nmf_h) / rowSums(t(nmf_h))
    ## Then with this - we define "cluster core" as set of samples with cluster membership > 0.5
    # clustercores <- nmf_h_membership_scores > 0.5
    
    ## Extract imporatnt features
    feature_scores <- featureScore(nmfout_object) # gives a score for each gene from ~0 - ~1 
    driving_features <- extractFeatures(nmfout_object) # for each rank, outputs the most important genes (pretty sure this is right)
    ## Get the features for the 3rd cluster
    # feature_scores1[driving_features[[3]]] #CD177 and MMP9 - pretty cool actually.......
    # rownames(nmfout)[driving_features[[4]]]
    
    ## Outlist
    return(list(nmfobject = nmfout_object,
                nmf_w_mat = nmf_w,
                nmf_h_mat = nmf_h,
                nmf_h_membership_scores = nmf_h_membership_scores,
                nmf_featurescores = feature_scores,
                nmf_drivingfeatures = driving_features))
    
}




## Cbind fill function
cbind.fill <- function(...){
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

