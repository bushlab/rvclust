#' @include clustering.r
NULL

#' Collapses features into composite feature using the OR operator
#'
#' Performs a standard collapsing test, which collapses the features in a bin
#' by taking the OR of all features. For each subject, the composite feature
#' is represented as positive if any of the individual events have occurred.
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with collapsed data (composite features)
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
#' @seealso \code{\link{pamk}}
collapsing.test <- function(rv) {

	variants   <- rv$variants
 	rv$composite.features <- collapse.clusters(variants,rv$clusters,rv$data$ped)
 	return(rv)
}