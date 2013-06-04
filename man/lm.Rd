\name{lm}
\alias{lm}
\title{Linear Regression across clusters for rvclust (wrapper for stats::lm)}
\usage{
  lm(rv, label.by = NA, min.fit = 0, phen = "PHENOTYPE")
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  rvclustobject with linear regression results
}
\description{
  ----------------------------------------------------------------------
  '#\cr This module wraps stats::lm to perform linear
  regression across \cr collapsed clusters. These wrappers
  are not strictly necessary, \cr but do provide automated
  support for cluster fitness thresholds, \cr covariates,
  and calculation of statistical effect and significance.
  \cr
  ----------------------------------------------------------------------
  '#\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}

