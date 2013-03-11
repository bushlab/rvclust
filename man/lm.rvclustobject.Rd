\name{lm.rvclustobject}
\alias{lm.rvclustobject}
\title{Linear Regression across clusters for rvclust (wrapper for stats::lm)}
\usage{
  \method{lm}{rvclustobject} (rv)
}
\arguments{
  \item{rvclustobject}{}
}
\value{
  rvclustobject with linear regression results
}
\description{
  Performs linear regression by wrapping the stats::lm\cr
  function. Using cluster assignemnt, PEDMAP data is
  collapsed\cr into a boolean variable indicating the
  presence or absence\cr of a variant in that cluster.\cr
  \cr If burden testing has been specified, the boolean is
  replaced\cr with the total number of variants in each
  cluster. \cr Only those clusters with a fitness exceeding
  the minimum fitness\cr threshold are tested. They should
  not be included in multiple\cr test correction. If no
  minimum fitness has been threshold was\cr was specified,
  all clusters are tested for association. \cr If any
  covariates were specified when creating the
  rvclustobject,\cr those covariates will be included in
  the association test. \cr Association results will be
  appended to the clusterinfo object.
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}

  \code{\link{pamk}}

  \code{\link{hpower}}
}

