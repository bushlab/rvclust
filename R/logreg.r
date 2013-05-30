#' @include clustering.r
NULL

#' Logistic Regression across clusters for rvclust (wrapper for rms::lrm)
#' 
#' ---------------------------------------------------------------------- '#\cr
#' This module wraps rms::lrm to perform logistic regression across \cr
#' collapsed clusters. These wrappers are not strictly necessary, \cr
#' but do provide automated support for cluster fitness thresholds, \cr
#' covariates, and calculation of statistical effect and significance. \cr
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with logistic regression results
logreg <- function(rv,min.fit=0.0,phen="PHENOTYPE") {
  rarevariants <- rv$variants
  clusterinfo <- rv$clusterinfo
  cov.dat <- rv$covariates
  phen.dat <- rv$observations
  collapsed.dat <- rv$collapsed
  
  # If the clustering algorithm recorded fitness, filter clusters below 
  # the minimum fitness threshold
  if ("FIT" %in% names(clusterinfo)) {
    clusterinfo <- clusterinfo[clusterinfo$FIT>min.fit,] }
  
  # Extract the unique cluster IDs
  clusterids <- as.character(unique(clusterinfo$CLUSTERID))  
    
  # Add the covariates
  covariates <- NA
  if (!any(is.na(cov.dat))) {
    covariates <- names(cov.dat)
    covariates <- covariates[3:length(covariates)]  # Filter out IDs
    collapsed.dat <- merge(collapsed.dat,cov.dat,by='FID')
  }
  phenotypes <- NA
  # Add additional phenotypes if provided
  if (!any(is.na(phen.dat))) {
    phenotypes <- names(phen.dat)
    phenotypes <- phenotypes[3:length(phenotypes)]  # Filter out IDs
    collapsed.dat <- merge(x=collapsed.dat,y=phen.dat,by='FID')
  }
  
  # Run a logistic regression for each CLUSTERID
  models <- sapply(clusterids,function(id){
    predictors <- paste("collapsed.dat$\"X",id,"\"",sep='')
    if (!any(is.na(covariates))) {
      predictors <- paste(append(predictors,covariates),collapse='+')}
    f <- as.formula(paste(phen," ~ ",predictors,sep=''))
    cluster.logreg <- rms::lrm(f, data=collapsed.dat)
    return(cluster.logreg)
  },simplify=FALSE)
  
  pvalues <- sapply(models,function(model) {model$stats["P"]})
  r.squareds <- sapply(models,function(model) {model$stats["R2"]})

  clusterinfo$PVALUE <- pvalues
  clusterinfo$R.SQUARED <- r.squareds

  rv$clusterinfo <- clusterinfo
  return(rv)
}
