#' Automated test for statistical association modules
#'
#' ---------------------------------------------------------------------- '#\cr
#' Use this function to test statistical association modules.
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param association function
#' @return validity
test.association <- function(ASSOCIATION.FUN) {
	data('test.rv')

	# Cluster with rvclust::pamk
	rv <- pamk(rv,label.by='rsid',cluster.by='pos')
	rv <- collapsing.test(rv)

	# Check for runtime errors
	rv <- tryCatch({
		ASSOCIATION.FUN(rv,label.by='FID')
		},error=function(err){
			stop(paste("There was a runtime error.",err,sep='\n'))
		})

	# Check for valid rvclustobject
	if (!(valid.rvclustobject(rv))) {
		stop("The returned rvclustobject is invalid.")
	}

	# Everything looks good
	writeLines("Association module is valid.")
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
	if (!("collapsed" %in% names(rv))) {
		missing.fields <- append(missing.fields,"collapsed")}
	if (length(missing.fields) > 0) {
		fields <- paste(missing.fields,collapse=',')
		writeLines(paste("Missing fields",fields,sep=': '))
		return(FALSE)
	}
	return(TRUE)
}
}