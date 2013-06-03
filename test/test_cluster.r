#' Automated test for clustering modules
#'
#' ---------------------------------------------------------------------- '#\cr
#' Use this function to test cluster modules.
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return annotated rvclustobject
test.cluster <- function(CLUSTER.FUN) {
	print("Functionality not yet implemented.")
	load('test.rv.rda')

	# Check for runtime errors
	rv <- tryCatch({
		CLUSTER.FUN(rv,label.by='rsid',cluster.by='pos')
		},error=function(err){
			stop(paste("There was an error during clustering.",err,sep='\n'))
		})

	# Check for valid rvclustobject
	if (!(valid.rvclustobject(rv))) {
		stop("The returned rvclustobject is invalid.")
	}

	# Everything looks good
	writeLines("Clustering module is valid.")
}

valid.rvclustobject <- function(rv) {
	missing.fields <- c()
	if (!("observations" %in% names(rv))) {
		missing.fields <- append(missing.fields,"observations")}
	if (!("variants" %in% names(rv))) {
		missing.fields <- append(missing.fields,"variants")}
	if (!("covariates" %in% names(rv))) {
		missing.fields <- append(missing.fields,"covariates")}
	if (!("outcomes" %in% names(rv))) {
		missing.fields <- append(missing.fields,"outcomes")}
	if (!("clusters" %in% names(rv))) {
		missing.fields <- append(missing.fields,"clusters")}
	if (!("clusterinfo" %in% names(rv))) {
		missing.fields <- append(missing.fields,"clusterinfo")}
	if (length(missing.fields) > 0) {
		fields <- paste(missing.fields,collapse=',')
		stop(paste("Missing fields",fields,sep=': '))
	}
}