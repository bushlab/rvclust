#' @include clustering.r
NULL

#' K-Medoids Clustering for RVCLUST (wrapper for fpc::pamk)
#'
#' ---------------------------------------------------------------------- '#\cr
#' This module wraps fpc::pamk to perform k-Medoids clustering. To apply \cr
#' constraints, specify a set of columns that must be cluster homogenous. \cr
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return clustered rvclustobject
pamk <- function(rv,label.by=NA,cluster.by=NA,constrain.by=NA) {
  variants <- rv$variants

  # If a constraints are specified, force
  # homogenous clusters wrt constrain.by
  if (!is.na(constrain.by)) {
    f <- variants[,constrain.by]
    variants$CLUSTERID <- rep(0,length(f))
    df.matrix <- split(variants,f)
  }
  else {
    df.matrix <- list(variants)
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
  
  # Rebuild the variants dataframe (unsplit) (1-op for unconstrained)
  variants <- df.matrix[[1]]
  if (length(df.matrix) > 1) {
    for (i in 2:length(df.matrix)) {
      variants <- rbind(variants,df.matrix[[i]])}
  }
  
  # Create a vector of clustered SNPs
  clusters         <- list()
  length(clusters) <- k
  clusters         <- sapply(clusters,function(x){rep('',0)})
  for (i in 1:nrow(variants)) {
    cluster.id <- variants[i,"CLUSTERID"]
    snp        <- as.character(variants[i,label.by])
    res        <- append(clusters[[cluster.id]],snp)
    clusters[[cluster.id]] <- res
  }

  # Create a cluster info dataframe
  clusterinfo    <- create.cluster.dat(clusters,variants)

  rv$clusters    <- clusters
  rv$clusterinfo <- clusterinfo
  rv$variants    <- variants
  return(rv)
}
