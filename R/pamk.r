#' @include clustering.r
NULL

#require(fpc)
#source('clustering.r')
#' Clusters rare variants using fpc::pamk
#'
#' Wrapper for fpc::pamk standardizing input and output\cr
#' for rvclust.
#'
#' @export
#' @param rarevariants rare variant data frame
#' @param raw.dat genotype/phenotype data frame
#' @param vars list of strings identifying the cluster columns
#' @param key string identifying the label column
#' @return cluster object containing\cr
#'  \tabular{ll}{
#'  snp.clusters: \tab list containing SNPs grouped by cluster\cr
#'  cluster.info: \tab data frame containing information on each cluster\cr
#'  rarevariants: \tab rare variant data frame with cluster assignment added
#'  }
#' @seealso \code{\link{init}}
#' @seealso \code{\link{annotate}}
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @examples
pamk <- function(rv) {
  rarevariants <- rv$variants
  raw.dat <- rv$data$ped
  
  # Split the data across chromatin state
  f <- rarevariants$STATE
  rarevariants$CLUSTERID <- rep(0,length(f))
  df.matrix <- split(rarevariants,f)
  
  # Cluster for each chromatin state
  k <- 0
  for (i in 1:length(df.matrix)) {
    df <- df.matrix[i][[1]]
    krange <- 2:min(15,length(df$POS)-1)
    clusters <- fpc::pamk(df$POS,krange=krange)
    df$CLUSTERID <- clusters$pamobject$clustering+k
    k <- clusters$nc+k
    df.matrix[i][[1]] <- df
  }
  
  # Rebuild the rarevariants dataframe (unsplit)
  rarevariants <- df.matrix[1][[1]]
  for (i in 2:length(df.matrix)) {
    rarevariants <- rbind(rarevariants,df.matrix[i][[1]])}
  
  # Create a cluster info dataframe
  cluster.info <- create.cluster.dat(k,rarevariants)
  
  # Create a vector of clustered SNPs
  snp.clusters <- list()
  length(snp.clusters) <- k
  snp.clusters <- sapply(snp.clusters,function(x){rep('',0)})
  for (i in 1:nrow(rarevariants)) {
    cluster.id <- rarevariants[i,"CLUSTERID"]
    snp <- as.character(rarevariants[i,"SNP"])
    res <- append(snp.clusters[[cluster.id]],snp)
    snp.clusters[[cluster.id]] <- res
  }
 
  rv$clusters <- snp.clusters
  rv$clusterinfo <- cluster.info
  rv$variants <- rarevariants
  return(rv)
}