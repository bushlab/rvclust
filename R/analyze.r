#' @include clustering.r
NULL

#' Analyzes clusters for statistical significance
#' 
#' Analyzes clusters for statistical significance using\cr
#' either a collapsing window or burden test. Testing only\cr
#' those clusters with a specified minimum fitness is encouraged,\cr
#' but not required.
#'
#' @export
#' @param rarevariants rare variant data frame
#' @param clusterinfo data frame containing information on each cluster
#' @param raw.dat genotype/phenotype data frame
#' @param cov.dat covariate data frame
#' @param burden boolean indicating whether to use burden testing
#' @param min.fit real specifying the minimum cluster fitness to evaluate
#' @return Cluster data frame with statistical results columns
#' @seealso \code{\link{init}}
#' @seealso \code{\link{annotate}}
#' @seealso \code{\link{pamk}}
#' @seealso \code{\link{hpower}}
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @examples
#' analyze(rarevariants,clusterinfo,raw.dat)
#' analyze(rarevariants,clusterinfo,raw.dat,cov.dat)
#' analyze(rarevariants,clusterinfo,raw.dat,burden=TRUE)
#' analyze(rarevariants,clusterinfo,raw.dat,min.fit=0.5)
#' analyze(rarevariants,clusterinfo,raw.dat,cov.dat,burden=TRUE,min.fit=0.5)
analyze <- function(rv,burden=FALSE,min.fit=0) {
  rarevariants <- rv$variants
  clusterinfo <- rv$clusterinfo
  raw.dat <- rv$data$ped
  cov.dat <- rv$data$cov
  
  # Reduce to those clusters meeting minimum fitness
  clusterinfo <- clusterinfo[clusterinfo$FIT>min.fit,]
  
  # Extract the unique cluster IDs
  clusterids <- as.character(unique(clusterinfo$CLUSTERID))
  
  # Collapse the raw dataset into clusters
  collapsed.dat <- collapse.clusters(rarevariants,raw.dat,burden)
  
  # Write the collapsed.dat to file for later reference
  #write.table(collapsed.dat,paste('temp/',data.fname,"_collapsed.raw",sep=''),col.names=TRUE,row.names=FALSE)
  
  # Add the covariates
  covariates <- NA
  if (!is.na(cov.dat)) {
    covariates <- names(cov.dat)
    collapsed.dat <- merge(collapsed.dat,cov.dat,by.x=FID,by.y=cov.dat[1])
  }
  
  # Run an anova for each CLUSTERID
  models <- sapply(clusterids,function(id){
    predictors <- paste("collapsed.dat$\"",id,"\"",sep='')
    if (!is.na(covariates)) {
      predictors <- paste(append(predictors,covariates),collapse='+')}
    f <- as.formula(paste("PHENOTYPE ~ ",predictors,sep=''))
    cluster.lm <- lm(f, data=collapsed.dat)
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