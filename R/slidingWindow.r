#' @include clustering.r
NULL

#' Sliding Window Analysis for rvclust (wrapper for fpc::pamk)
#'
#' Performs sliding window analysis. Included for comparison.\cr
#' \cr
#' Cluster assignment will be appended to the variants object.\cr
#' The clusters object will be populated with SNP lists for each\cr
#' window. The clusterinfo object will be populated with the\cr
#' following information about each window:\cr
#' 1. Chromosome\cr
#' 2. Minimum Base Pair\cr
#' 3. Maximum Base Pair\cr
#' 4. Size\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return clustered rvclustobject
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
sliding.window <- function(rv,win.size=5000,label.by=NA,cluster.by=NA,constrain.by=NA) {
	variants <- rv$variants
	range <- win.size/2

	# Create a vector of clustered SNPs
	snp.clusters <- sapply(variants[,cluster.by],function(x,variants,range)
		{variants[abs(x-variants[,cluster.by])<range,][,label.by]},
		variants=variants,range=range)

	# Create a cluster info dataframe
	cluster.info <- create.cluster.dat(snp.clusters,variants)

	rv$clusters <- snp.clusters
 	rv$clusterinfo <- cluster.info
 	return(rv)
 }