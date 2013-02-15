\docType{package}
\name{rvclust-package}
\alias{rvclust}
\alias{rvclust-package}
\title{rvclust: Methods for clustering rare variants for statistical analysis}
\description{
  Clusters rare variants for statistical analysis
}
\details{
  \tabular{ll}{ Package: \tab rvclust\cr Type: \tab
  Package\cr Version: \tab 2.0\cr Date: \tab 2013-02-01\cr
  License: \tab GPL (>= 2)\cr LazyLoad: \tab yes\cr }

  A package for applying various clustering\cr algorithms
  to rare variant analysis with\cr either burden or
  collapsed testing.\cr \cr Suggested Workflow:\cr 1)
  Create an rvclustobject \code{\link{rvclustobject}}\cr 2)
  Annotate data with \code{\link{annotate}}\cr 3) Cluster
  data with \code{\link{pamk}} or \code{\link{hpower}}\cr
  4) Analyze clusters with
  \code{\link{analyze.clusters}}\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate.rvclustobject}}

  \code{\link{pamk.rvclustobject}}

  \code{\link{rvcluster.rvclustobject}}

  \code{\link{lm.rvclustobject}}
}
\keyword{package}

