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
#' @method sliding.window rvclustobject
#' @param rv rvclustobject
#' @return clustered rvclustobject
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
sliding.window <- function(rv,label.by=NA,cluster.by=NA,constrain.by=NA) {
	variants <- rv$variants
	range <- 2500

	# Create a vector of clustered SNPs
	snp.clusters <- sapply(variants$POS,function(x,variants,range)
		{variants[abs(x-variants$POS)<range,]$SNP},
		variants=variants,range=range)

	# Create a cluster info dataframe
	cluster.info <- create.cluster.dat(length(snp.clusters),variants)

	# Initialize empty lists to store cluster IDs
	rv$variants$CLUSTERID = rep(c(),nrow(rv$variants))

	# Given a snp and clusterid, update the rvclustobject
	identify <- function(snp,clusterid,rv) {
		clusters <- rv$variants[rv$variants$SNP==snp,]$CLUSTERID
		rv$variants[rv$variants$SNP==snp,]$CLUSTERID <- append(clusters,clusterid)
	}

	# Assign all cluster IDs
	for i in 1:length(snp.clusters) {
		sapply(snp.clusters[[i]],identify,clusterid=i,rv=rv)}

	rv$clusters <- snp.clusters
 	rv$clusterinfo <- cluster.info
 	return(rv)
 }