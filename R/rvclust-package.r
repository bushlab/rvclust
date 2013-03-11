#' Clusters rare variants for statistical analysis
#'
#' \tabular{ll}{
#' Package: \tab rvclust\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2013-02-01\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' rvclust provides an object oriented approach to rare variant clustering\cr
#' for the purposes of testing genetic association through collapsed window\cr
#' or burden testing.
#' \cr
#' Initialize an rvclustobject with PEDMAP file locations, covariates,\cr
#' and configuration settings. Annotate your data with GWAR using the\cr
#' annotate function. Choose a clustering method interface or write an\cr
#' interface to one of your own. Test for associations using one of the\cr
#' provided interfaces, or again, write your own.
#'
#'
#' @name rvclust-package
#' @aliases rvclust
#' @docType package
#' @title rvclust: Methods for clustering rare variants for statistical analysis
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @keywords package
#' @seealso \code{\link{rvclustobject}}
#' @seealso \code{\link{annotate.rvclustobject}}
#' @seealso \code{\link{pamk.rvclustobject}}
#' @seealso \code{\link{rvcluster.rvclustobject}}
#' @seealso \code{\link{lm.rvclustobject}}
#' @references Sivley, R. Michael, Fish, Alexandra E., Bush, William S. (2013). 
#'  Knowledge-constrained K-medoids Clustering of  Regulatory Rare
#'  Alleles for Burden Tests. Lecture Notes in Computer Science, 
#'  vol. 7833 (In Press)
NULL