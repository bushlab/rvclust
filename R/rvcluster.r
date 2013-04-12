#' @include clustering.r
NULL

#' Hierarchical Clustering using Data Entropy for rvclust
#'
#' Top down hierarchical clustering where each level is partitioned\cr
#' into two clusters using k-Medoids (cluster::pamk) using base\cr
#' position and any annotations provided by the annotate function.\cr
#' \cr
#' The frontier is determined using data entropy as a proxy for\cr
#' statistical power. If a cluster has higher data entropy than both\cr
#' of its child clusters, the children are discarded and the parent\cr
#' is chosen as the leaf node.\cr
#' \cr
#' This strategy results in a frontier of two types of clusters.\cr
#' 1. Fit clusters chosen for their high data entropy\cr
#' 2. Unfit clusters removed from a parent to create a more fit sibling\cr
#' It is important when using this method to set a minimum fitness\cr
#' threshold so that the unfit clusters are not included in the\cr
#' analysis.\cr
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
#' 6. Annotations (multi-column)
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @method rvcluster rvclustobject
#' @param rv rvclustobject
#' @return clustered rvclustobject
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
rvcluster <- function(rv) {
  rarevariants <- rv$variants
  raw.dat <- rv$data$ped
  key <- "SNP"
  vars <- c("POS",rv$annotations)
  
  rarevariants$CLUSTERID <- rep(0,length(rarevariants$POS))
  
  vars <- append(vars,key)
  cluster.tree <- entropyClust(rarevariants[,vars],key,raw.dat,'entropy')

  # Function to parse the tree and return the frontier as a list of lists
  get.frontier <- function(node) {
    if (is.na(node[[2]]) | is.na(node[[3]])) {
      return(list(node[[1]]))}
    else {
      left  <- get.frontier(node[[2]])
      right <- get.frontier(node[[3]])
      return(c(left, right))
    }
  }
  frontier <- get.frontier(cluster.tree)

  # Function to parse the frontier and return the fitnesses as a list
  get.fitness <- function(node) {
    if (is.na(node[[2]]) | is.na(node[[3]])) {
      return(list(node[[4]]))}
    else {
      left  <- get.fitness(node[[2]])
      right <- get.fitness(node[[3]])
      return(c(left, right))
    }
  }
  fitness <- unlist(get.fitness(cluster.tree))

  # Assign cluster IDs to each data frame
  assign.cluster <- function(frontier,rarevariants) {
    for (i in 1:length(frontier)) {
      for (snp in frontier[[i]]) {
        rarevariants[rarevariants$SNP==snp,]$CLUSTERID <- i
      }
    }
    return(rarevariants)
  }
  
  rarevariants <- assign.cluster(frontier,rarevariants)

  # Create a cluster info dataframe
  k <- length(frontier)
  cluster.info <- create.cluster.dat(k,rarevariants,NA,fitness,vars)
  
  rv$clusters <- frontier
  rv$clusterinfo <- cluster.info
  rv$variants <- rarevariants
  return(rv)
}
