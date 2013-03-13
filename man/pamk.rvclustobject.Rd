\name{pamk.rvclustobject}
\alias{pamk.rvclustobject}
\title{K-Medoids Clustering for rvclust (wrapper for fpc::pamk)}
\usage{
  \method{pamk}{rvclustobject} (rv)
}
\arguments{
  \item{rvclustobject}{}
}
\value{
  clustered rvclustobject
}
\description{
  Performs K-Medoids clustering by wrapping the fpc::pamk
  function.\cr \cr Cluster assignment will be appended to
  the data$ped object.\cr The clusters object will be
  populated with SNP lists for each\cr cluster. The
  clusterinfo object will be populated with the\cr
  following information about each cluster:\cr 1.
  Chromosome\cr 2. Minimum Base Pair\cr 3. Maximum Base
  Pair\cr 4. Size\cr 5. Fitness\cr 6. Annotations
  (multi-column)\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}
}

