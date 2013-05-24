#' @include clustering.r
NULL

#' Collapses features into composite feature using the OR operator
#'
#' Performs a standard burden test, which collapses the features in a bin
#' by taking the SUM of all features. For each subject, the composite
#' feature is represented by the burden of events occurring in the bin.
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with collapsed data (composite features)
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
#' @seealso \code{\link{pamk}}
burden.test <- function(rv) {

  variants   <- rv$variants
  rv$collapsed <- collapse.clusters(variants,rv$clusters,rv$observations,burden=TRUE)
  return(rv)
}