\name{rvcluster.rvclustobject}
\alias{rvcluster.rvclustobject}
\title{Hierarchical Clustering using Data Entropy for rvclust}
\usage{
  \method{rvcluster}{rvclustobject} (rv)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  clustered rvclustobject
}
\description{
  Top down hierarchical clustering where each level is
  partitioned\cr into two clusters using k-Medoids
  (cluster::pamk) using base\cr position and any
  annotations provided by the annotate function.\cr \cr The
  frontier is determined using data entropy as a proxy
  for\cr statistical power. If a cluster has higher data
  entropy than both\cr of its child clusters, the children
  are discarded and the parent\cr is chosen as the leaf
  node.\cr \cr This strategy results in a frontier of two
  types of clusters.\cr 1. Fit clusters chosen for their
  high data entropy\cr 2. Unfit clusters removed from a
  parent to create a more fit sibling\cr It is important
  when using this method to set a minimum fitness\cr
  threshold so that the unfit clusters are not included in
  the\cr analysis.\cr \cr Cluster assignment will be
  appended to the data$ped object.\cr The clusters object
  will be populated with SNP lists for each\cr cluster. The
  clusterinfo object will be populated with the\cr
  following information about each cluster:\cr 1.
  Chromosome\cr 2. Minimum Base Pair\cr 3. Maximum Base
  Pair\cr 4. Size\cr 5. Fitness\cr 6. Annotations
  (multi-column)
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}
}

