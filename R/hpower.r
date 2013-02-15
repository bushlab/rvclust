#' @include clustering.r
NULL

#' Clusters rare variants using rvclust custom hierarchical method
#'
#' Top down hierarchical cluster where each level is partitioned\cr
#' into two clusters using fpc::pamk across the specified columns.\cr
#' The fitness function, which decides whether to expand a cluster\cr
#' is a statistical power proxy provided by Dr. Eli Stahl, Emory University. 
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
#'  rarevariants: \tab rare variant data frame with cluster assignment added\cr
#'  }
#' @seealso \code{\link{init}}
#' @seealso \code{\link{annotate}}
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @examples
hpower <- function(rv,vars=c("POS","STATE"),key="SNP") {
  rarevariants <- rv$variants
  raw.dat <- rv$data$ped

  # Split the data across chromatin state
  #f <- rarevariants$STATE
  #rarevariants$CLUSTERID <- rep(0,length(f))
  #df.matrix <- split(rarevariants,f)
  
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
  cluster.info <- create.cluster.dat(k,rarevariants,NA,fitness)
  
  rv$clusters <- frontier
  rv$clusterinfo <- cluster.info
  rv$variants <- rarevariants
  return(rv)
}
