\name{sliding.window.rvclustobject}
\alias{sliding.window.rvclustobject}
\title{Sliding Window Analysis for rvclust (wrapper for fpc::pamk)}
\usage{
  \method{sliding.window}{rvclustobject} (rv)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  clustered rvclustobject
}
\description{
  Performs sliding window analysis. Included for
  comparison.\cr \cr Cluster assignment will be appended to
  the variants object.\cr The clusters object will be
  populated with SNP lists for each\cr window. The
  clusterinfo object will be populated with the\cr
  following information about each window:\cr 1.
  Chromosome\cr 2. Minimum Base Pair\cr 3. Maximum Base
  Pair\cr 4. Size\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}
}

