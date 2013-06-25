#' Performs single-marker linear regression on variants
#'
#' See above
#' 
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with variant-level linear regression results
single.marker.lm <- function(rv,label.by=NA,phen="PHENOTYPE") {
  variants    <- rv$variants
  cov.dat     <- rv$covariates
  phen.dat    <- rv$outcomes
  observations<- rv$observations

  # Add the covariates if provided
  covariates <- NA
  if (!any(is.na(cov.dat))) {
    covariates <- names(cov.dat)
    covariates <- covariates[3:length(covariates)]  # Filter out IDs
    observations <- merge(x=observations,y=cov.dat,by=label.by)
  }

  phenotypes <- NA
  # Add additional phenotypes if provided
  if (!any(is.na(phen.dat))) {
    observations <- merge(x=observations,y=phen.dat,by=label.by)
  }

  #FIXME Assumption made that the variants labels are in the first column
  variant.names <- variants[,1]
  models <- sapply(variant.names,function(variant){
    predictors <- variant
    if (!any(is.na(covariates))) {
      predictors <- paste(append(predictors,covariates),collapse='+')}
    f <- as.formula(paste(phen," ~ ",predictors,sep=''))
    variant.lm <- stats::lm(f, data=observations)
    variant.lm
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

  variants$PVALUE <- pvalues
  variants$EFFECT <- effects
  variants$R.SQUARED <- r.squareds
  variants$ADJ.R.SQUARED <- adj.r.squareds

  rv$variants <- variants
  return(rv)
}

#rvobj.test <- single.marker.lm(rvobj,label.by='FID')