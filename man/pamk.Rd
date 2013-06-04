\name{pamk}
\alias{pamk}
\title{K-Medoids Clustering for RVCLUST (wrapper for fpc::pamk)}
\usage{
  pamk(rv, label.by = NA, cluster.by = NA,
    constrain.by = NA)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  clustered rvclustobject
}
\description{
  ----------------------------------------------------------------------
  '#\cr This module wraps fpc::pamk to perform k-Medoids
  clustering. To apply \cr constraints, specify a set of
  columns that must be cluster homogenous. \cr
  ----------------------------------------------------------------------
  '#\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}

