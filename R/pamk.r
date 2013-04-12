#' @include clustering.r
NULL

#' K-Medoids Clustering for rvclust (wrapper for fpc::pamk)
#'
#' Performs K-Medoids clustering by wrapping the fpc::pamk function.\cr
#' \cr
#' Cluster assignment will be appended to the data$ped object.\cr
#' The clusters object will be populated with SNP lists for each\cr
#' cluster. The clusterinfo object will be populated with the\cr
#' following information about each cluster:\cr
#' 1. Chromosome\cr
#' 2. Minimum Base Pair\cr
#' 3. Maximum Base Pair\cr
#' 4. Size\cr
#' 5. Fitness\cr
#' 6. Annotations (multi-column)\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @method pamk rvclustobject
#' @param rv rvclustobject
#' @return clustered rvclustobject
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
pamk <- function(rv,label.by=NA,cluster.by=NA,constrain.by=NA) {
  rarevariants <- rv$variants
  raw.dat <- rv$data$ped
  
  # Split the data across chromatin state
  f <- rarevariants$CHROMATIN
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
