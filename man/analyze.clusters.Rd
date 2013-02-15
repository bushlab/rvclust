\name{analyze.clusters}
\alias{analyze.clusters}
\title{Analyzes clusters for statistical significance}
\usage{
  analyze.clusters(rv.dat, cluster.dat, raw.dat,
    cov.dat = NA, burden = FALSE, min.fit = 0)
}
\arguments{
  \item{rv.dat}{rare variant data frame}

  \item{cluster.dat}{data frame containing information on
  each cluster}

  \item{raw.dat}{genotype/phenotype data frame}

  \item{cov.dat}{covariate data frame}

  \item{burden}{boolean indicating whether to use burden
  testing}

  \item{min.fit}{real specifying the minimum cluster
  fitness to evaluate}
}
\value{
  Cluster data frame with statistical results columns
}
\description{
  Analyzes clusters for statistical significance using\cr
  either a collapsing window or burden test. Testing
  only\cr those clusters with a specified minimum fitness
  is encouraged,\cr but not required.
}
\examples{
analyze.clusters(rv.dat,cluster.dat,raw.dat)
analyze.clusters(rv.dat,cluster.dat,raw.dat,cov.dat)
analyze.clusters(rv.dat,cluster.dat,raw.dat,burden=TRUE)
analyze.clusters(rv.dat,cluster.dat,raw.dat,min.fit=0.5)
analyze.clusters(rv.dat,cluster.dat,raw.dat,cov.dat,burden=TRUE,min.fit=0.5)
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{init}}

  \code{\link{annotate}}

  \code{\link{pamk}}

  \code{\link{hpower}}
}

