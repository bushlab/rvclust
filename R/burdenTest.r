#' @include clustering.r
NULL

#' Collapses features into composite feature using the OR operator
#'
#' ---------------------------------------------------------------------- '#\cr
#' Performs a standard burden test, which collapses the features in a bin \cr
#' by taking the SUM of all features. For each subject, the composite \cr
#' feature is represented by the burden of events occurring in the bin. \cr
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return rvclustobject with collapsed data (composite features)
burden.test <- function(rv) {
  variants   <- rv$variants
  rv$collapsed <- collapse.clusters(variants,rv$clusters,rv$observations,burden=TRUE)
  return(rv)
}