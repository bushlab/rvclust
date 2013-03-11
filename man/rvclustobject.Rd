\name{rvclustobject}
\alias{rvclustobject}
\title{Constructor for rvclustobject}
\usage{
  rvclustobject(pedmap.path, pedmap.fname, cov.path = NA,
    cov.fname = NA, annotations = NA, burden = FALSE,
    min.fit = 0)
}
\arguments{
  \item{pedmap.path}{path to pedmap files}

  \item{pedmap.fname}{basename for pedmap files (no
  extension)}

  \item{cov.path}{path to covariates file [optional]}

  \item{cov.fname}{filename for covariates file (with
  extension) [optional]}

  \item{annotations}{list of desired annotations [default
  ALL]}

  \item{burden}{boolean indicating whether to burden test
  [default FALSE]}

  \item{min.fit}{minimum cluster fitness to test (0..1)
  [default 0.0]}
}
\value{
  rvclustobject containing\cr \tabular{ll}{ data: \tab
  ped/map/covariate data\cr variants: \tab variant
  data.frame\cr clusters: \tab empty list of clusters\cr
  clusterinfo: \tab empty list of cluster info\cr }
}
\description{
  Constructor for rvclustobject. Loads ped/map data from
  the\cr specified path and filename. If a covariate path
  and\cr filename are specified, the covariates will be
  loaded as well.\cr Burden testing can be enabled by
  setting the burden flag to true.\cr A threshold can also
  be set such that only clusters with fitness\cr exceeding
  the threshold are tested for association.\cr \cr The
  rvclust package contains wrappers for various
  clustering\cr methods and statistical tests. See package
  documentation for more\cr information on these wrappers.
  Simply pass the rvclustobject to any\cr method with an
  rvclust wrapper and an updated rvclustobject will be
  returned.
}
\note{
  Data must be in PEDMAP format

  Covariates - One column for subject; One column per
  covariate
}
\examples{
rvclustobject('/path/to','pedmap')
rvclustobject(...,cov.path='/path/to',cov.fname='covariates.csv')
rvclustobject(...,annotations=c("CHROMATIN","OTHER"))
rvclustobject(...,burden=TRUE)
rvclustobject(...,min.fit=0.8)
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}

