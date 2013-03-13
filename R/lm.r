#' @include clustering.r
NULL

#' lm generic
#'
#' @export
lm <- function(x) {UseMethod("lm",x)}
#' Linear Regression across clusters for rvclust (wrapper for stats::lm)
#' 
#' Performs linear regression by wrapping the stats::lm\cr
#' function. Using cluster assignemnt, PEDMAP data is collapsed\cr
#' into a boolean variable indicating the presence or absence\cr
#' of a variant in that cluster.\cr
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
#' @method lm rvclustobject
#' @param rvclustobject
#' @return rvclustobject with linear regression results
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
#' @seealso \code{\link{pamk}}
#' @seealso \code{\link{hpower}}
lm.rvclustobject <- function(rv) {
  rarevariants <- rv$variants
  clusterinfo <- rv$clusterinfo
  raw.dat <- rv$data$ped
  cov.dat <- rv$data$cov
  burden <- rv$burden
  min.fit <- rv$min.fit
  
  # Reduce to those clusters meeting minimum fitness if
  # the clustering algorithm recorded fitness
  if ("FIT" %in% names(clusterinfo)) {
    clusterinfo <- clusterinfo[clusterinfo$FIT>min.fit,] }
  
  # Extract the unique cluster IDs
  clusterids <- as.character(unique(clusterinfo$CLUSTERID))
  
  # Collapse the raw dataset into clusters
  collapsed.dat <- collapse.clusters(rarevariants,raw.dat,burden)
    
  # Add the covariates
  covariates <- NA
  if (!is.na(cov.dat)) {
    covariates <- names(cov.dat)
    collapsed.dat <- merge(collapsed.dat,cov.dat,by.x=FID,by.y=cov.dat[1])
  }
  
  # Run a linear regression for each CLUSTERID
  models <- sapply(clusterids,function(id){
    predictors <- paste("collapsed.dat$\"",id,"\"",sep='')
    if (!is.na(covariates)) {
      predictors <- paste(append(predictors,covariates),collapse='+')}
    f <- as.formula(paste("PHENOTYPE ~ ",predictors,sep=''))
    cluster.lm <- stats::lm(f, data=collapsed.dat)
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