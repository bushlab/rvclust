\name{logreg}
\alias{logreg}
\title{Logistic Regression across clusters for rvclust (wrapper for rms::lrm)}
\usage{
  logreg(rv, label.by = NA, min.fit = 0,
    phen = "PHENOTYPE")
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  rvclustobject with logistic regression results
}
\description{
  ----------------------------------------------------------------------
  '#\cr This module wraps rms::lrm to perform logistic
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

