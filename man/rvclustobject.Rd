\name{rvclustobject}
\alias{rvclustobject}
\title{Constructor for rvclustobject}
\usage{
  rvclustobject(observation.file = NA, variant.file = NA,
    cov.file = NA, outcome.file = NA, max.freq = 100)
}
\arguments{
  \item{observation.file}{observation data (or PED)}

  \item{variant.file}{variant data (or MAP)}

  \item{cov.file}{covariate data}

  \item{outcome.file}{outcome data}
}
\value{
  rvclustobject with fields\cr \tabular{ll}{ variants\cr
  observations\cr covariates\cr outcomes\cr clusters\cr
  clusterinfo\cr collapsed\cr }
}
\description{
  ----------------------------------------------------------------------
  '#\cr RVCLUST functionality is modularized for easy
  extension development. \cr Each module is standardized to
  accept an rvclustobject, which is \cr initialized with
  variant, observation, covariate, and outcome data \cr
  files. \cr \cr Variant annotations may be applied by
  \code{\link{annotate}}. \cr \cr RVCLUST was originally
  designed for the analysis of rare genetic \cr variants.
  In line with its original intent, PEDMAP files may be \cr
  used as input as an alternative to the standardized input
  format. \cr This functionality is dependent on PLINK. \cr
  ----------------------------------------------------------------------
  '#\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}

