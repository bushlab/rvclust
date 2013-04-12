#' @include clustering.r
NULL

#' burden.test generic
#'
#' @export
#' @param rv rvclustobject
burden.test <- function(rv) {UseMethod("burden.test",rv)}
#' Collapses features into composite feature using the OR operator
#'
#' Performs a standard burden test, which collapses the features in a bin
#' by taking the SUM of all features. For each subject, the composite
#' feature is represented by the burden of events occurring in the bin.
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @method burden.test rvclustobject
#' @param rv rvclustobject
#' @return rvclustobject with collapsed data (composite features)
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate}}
#' @seealso \code{\link{pamk}}
burden.test.rvclustobject <- function(rv) {

  variants   <- rv$variants
  rv$composite.features <- collapse.clusters(variants,rv$data$ped,burden=TRUE)
  return(rv)
}