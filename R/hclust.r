#' @include clustering.r
NULL

#' Clusters rare variants using stats::hclust
#'
#' Wrapper for stats::hclust standardizing input and output\cr
#' for rvclust. Frontier is computed using a statistical power
#' proxy. User may specifiy an alternative objective function.
#' If an objective function is specified, it should accept
#' a vector of values as its first parameter. Any arguments
#' to rvclust::hclust following the objective function
#' will be passed as additional parameters to the function.
#'
#' @param rv.dat rare variant data frame
#' @param raw.dat genotype/phenotype data frame
#' @param vars list of strings identifying the cluster columns
#' @param key string identifying the label column
#' @return cluster object containing\cr
#'  \tabular{ll}{
#'  snp.clusters: \tab list containing SNPs grouped by cluster\cr
#'  cluster.info: \tab data frame containing information on each cluster\cr
#'  rv.dat: \tab rare variant data frame with cluster assignment added
#'  }
#' @seealso \code{\link{init}}
#' @seealso \code{\link{annotate}}
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
hclust <- function(rv.dat,raw.dat,vars=c("POS","STATE"),key="SNP",objfun=power,...) {
  
  diss <- dist(rv.dat[,vars])
  cluster.tree <- stats::hclust(diss,krange=krange)
 
  return(list("snp.clusters"=snp.clusters,"cluster.info"=cluster.info,"rv.dat"=rv.dat))
}