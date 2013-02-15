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

  rvclust provides an object oriented approach to rare
  variant clustering\cr for the purposes of testing genetic
  association through collapsed window\cr or burden
  testing. \cr Initialize an rvclustobject with PEDMAP file
  locations, covariates,\cr and configuration settings.
  Annotate your data with GWAR using the\cr annotate
  function. Choose a clustering method interface or write
  an\cr interface to one of your own. Test for associations
  using one of the\cr provided interfaces, or again, write
  your own.
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\references{
  Sivley, R. Michael, Fish, Alexandra E., Bush, William S.
  (2013).  Knowledge-constrained K-medoids Clustering of
  Regulatory Rare Alleles for Burden Tests. Lecture Notes
  in Computer Science, vol. 7833 (In Press)
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate.rvclustobject}}

  \code{\link{pamk.rvclustobject}}

  \code{\link{rvcluster.rvclustobject}}

  \code{\link{lm.rvclustobject}}
}
\keyword{package}

