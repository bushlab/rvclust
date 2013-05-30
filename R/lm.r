#' Linear Regression across clusters for rvclust (wrapper for stats::lm)
#' 
#' ---------------------------------------------------------------------- '#\cr
#' This module wraps stats::lm to perform linear regression across \cr
#' collapsed clusters. These wrappers are not strictly necessary, \cr
#' but do provide automated support for cluster fitness thresholds, \cr
#' covariates, and calculation of statistical effect and significance. \cr
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with linear regression results
lm <- function(rv,label.by=NA,min.fit=0.0,phen="PHENOTYPE") {
  clusterinfo <- rv$clusterinfo
  cov.dat     <- rv$covariates
  phen.dat    <- rv$outcomes
  collapsed   <- rv$collapsed
  
  # If the clustering algorithm recorded fitness, filter clusters below 
  # the minimum fitness threshold
  if ("FIT" %in% names(clusterinfo)) {
    clusterinfo <- clusterinfo[clusterinfo$FIT>min.fit,] }
  
  # Extract the unique cluster IDs
  clusterids <- as.character(unique(clusterinfo$CLUSTERID))  
    
  # Add the covariates if provided
  covariates <- NA
  if (!any(is.na(cov.dat))) {
    covariates <- names(cov.dat)
    covariates <- covariates[3:length(covariates)]  # Filter out IDs
    collapsed <- merge(x=collapsed,y=cov.dat,by=label.by)
  }
  phenotypes <- NA
  # Add additional phenotypes if provided
  if (!any(is.na(phen.dat))) {
    phenotypes <- names(phen.dat)
    phenotypes <- phenotypes[3:length(phenotypes)]  # Filter out IDs
    collapsed <- merge(x=collapsed,y=phen.dat,by=label.by)
  }
  
  # Run a linear regression for each CLUSTERID
  models <- sapply(clusterids,function(id){
    predictors <- paste("collapsed$\"X",id,"\"",sep='')
    if (!any(is.na(covariates))) {
      predictors <- paste(append(predictors,covariates),collapse='+')}
    f <- as.formula(paste(phen," ~ ",predictors,sep=''))
    cluster.lm <- stats::lm(f, data=collapsed)
    cluster.lm
  },simplify=FALSE)
  
  pvalues <- sapply(models,function(model) {anova(model)$"Pr(>F)"[1]})
  effects <- sapply(models,function(model) {
    my.anova <- anova(model)
    ss  <- my.anova$"Sum Sq"
    pes <- ss/(ss+ss[length(ss)])
    pes <- pes[1]
    pes
  })
  r.squareds <- sapply(models,function(model) {summary(model)$"r.squared"})
  adj.r.squareds <- sapply(models,function(model) {summary(model)$"adj.r.squared"})

  clusterinfo$PVALUE <- pvalues
  clusterinfo$EFFECT <- effects
  clusterinfo$R.SQUARED <- r.squareds
  clusterinfo$ADJ.R.SQUARED <- adj.r.squareds

  rv$clusterinfo <- clusterinfo
  return(rv)
}
