#' @include clustering.r
NULL

#' Logistic Regression across clusters for rvclust (wrapper for rms::lrm)
#' 
#' Performs logistic regression by wrapping the rms::lrm\cr
#' function with family="binomial". Using cluster assignment,\cr
#' PEDMAP data is collapsed into a boolean variable indicating\cr
#' the presence or absence of a variant in that cluster.\cr
#' \cr
#' If burden testing has been specified, the boolean is replaced\cr
#' with the total number of variants in each cluster.
#' \cr
#' Only those clusters with a fitness exceeding the minimum fitness\cr
#' threshold are tested. They should not be included in multiple\cr
#' test correction. If no minimum fitness has been threshold was\cr
#' was specified, all clusters are tested for association.
#' \cr
#' If any covariates were specified when creating the rvclustobject,\cr
#' those covariates will be included in the association test.
#' \cr
#' Association results will be appended to the clusterinfo object.
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with logistic regression results
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
#' @seealso \code{\link{pamk}}
logreg <- function(rv,min.fit=0.0) {
  rarevariants <- rv$variants
  clusterinfo <- rv$clusterinfo
  cov.dat <- rv$data$cov
  
  # If the clustering algorithm recorded fitness:
  # Reduce to those clusters meeting minimum fitness
  if ("FIT" %in% names(clusterinfo)) {
    clusterinfo <- clusterinfo[clusterinfo$FIT>min.fit,] }
  
  # Extract the unique cluster IDs
  clusterids <- as.character(unique(clusterinfo$CLUSTERID))
  
  # Pull the clustered data
  collapsed.dat <- rv$composite.features
    
  # Add the covariates
  covariates <- NA
  if (!is.na(cov.dat)) {
    covariates <- names(cov.dat)
    collapsed.dat <- merge(collapsed.dat,cov.dat,by.x=FID,by.y=cov.dat[1])
  }
  
  # Run a logistic regression for each CLUSTERID
  models <- sapply(clusterids,function(id){
    predictors <- paste("collapsed.dat$\"X",id,"\"",sep='')
    if (!is.na(covariates)) {
      predictors <- paste(append(predictors,covariates),collapse='+')}
    f <- as.formula(paste("PHENOTYPE ~ ",predictors,sep=''))
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
