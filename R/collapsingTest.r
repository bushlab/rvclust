#' @include clustering.r
NULL

#' Collapses features into composite feature using the OR operator
#'
#' ---------------------------------------------------------------------- '#\cr
#' Performs a standard collapsing test, which collapses the features in a \cr
#' bin by taking the OR of all features. For each subject, the composite \cr
#' feature is represented as positive if any of the individual events \cr
#' have occurred. \cr
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with collapsed data
collapsing.test <- function(rv) {
	variants   <- rv$variants
 	rv$collapsed <- collapse.clusters(variants,rv$clusters,rv$observations)
 	return(rv)
}