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
#' @param rv rvclustobject
#' @return clustered rvclustobject
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
pamk <- function(rv,label.by=NA,cluster.by=NA,constrain.by=NA) {
  rarevariants <- rv$variants
  raw.dat <- rv$data$ped
  
  if (is.na(label.by)) {label.by <- "SNP"}
  if (is.na(cluster.by)) {cluster.by <- "POS"}

  # If a constraint variable is specified, force
  # homogenous clusters wrt constrain.by
  if (!is.na(constrain.by)) {
    f <- rarevariants[,constrain.by]
    rarevariants$CLUSTERID <- rep(0,length(f))
    df.matrix <- split(rarevariants,f)
  }
  else {
    df.matrix <- list(rarevariants)
  }
  
  # Cluster for each constrained data frame
  k <- 0
  for (i in 1:length(df.matrix)) {
    df <- df.matrix[[i]]
    krange <- 2:min(15,length(df$POS)-1)
    clusters <- fpc::pamk(df[,cluster.by],krange=krange)
    df$CLUSTERID <- clusters$pamobject$clustering+k
    k <- clusters$nc+k
    df.matrix[[i]] <- df
  }
  
  # Rebuild the rarevariants dataframe (unsplit)
  rarevariants <- df.matrix[[1]]
  if (length(df.matrix) > 1) {
    for (i in 2:length(df.matrix)) {
      rarevariants <- rbind(rarevariants,df.matrix[[i]])}
  }
  
  # Create a vector of clustered SNPs
  snp.clusters <- list()
  length(snp.clusters) <- k
  snp.clusters <- sapply(snp.clusters,function(x){rep('',0)})
  for (i in 1:nrow(rarevariants)) {
    cluster.id <- rarevariants[i,"CLUSTERID"]
    snp <- as.character(rarevariants[i,label.by])
    res <- append(snp.clusters[[cluster.id]],snp)
    snp.clusters[[cluster.id]] <- res
  }

  # Create a cluster info dataframe
  cluster.info <- create.cluster.dat(snp.clusters,rarevariants)

  rv$clusters <- snp.clusters
  rv$clusterinfo <- cluster.info
  rv$variants <- rarevariants
  return(rv)
}
