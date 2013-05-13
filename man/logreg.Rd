\name{logreg}
\alias{logreg}
\title{Logistic Regression across clusters for rvclust (wrapper for rms::lrm)}
\usage{
  logreg(rv, min.fit = 0)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  rvclustobject with logistic regression results
}
\description{
  Performs logistic regression by wrapping the rms::lrm\cr
  function with family="binomial". Using cluster
  assignment,\cr PEDMAP data is collapsed into a boolean
  variable indicating\cr the presence or absence of a
  variant in that cluster.\cr \cr If burden testing has
  been specified, the boolean is replaced\cr with the total
  number of variants in each cluster. \cr Only those
  clusters with a fitness exceeding the minimum fitness\cr
  threshold are tested. They should not be included in
  multiple\cr test correction. If no minimum fitness has
  been threshold was\cr was specified, all clusters are
  tested for association. \cr If any covariates were
  specified when creating the rvclustobject,\cr those
  covariates will be included in the association test. \cr
  Association results will be appended to the clusterinfo
  object.
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}

  \code{\link{pamk}}
}

