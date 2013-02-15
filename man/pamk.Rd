\name{pamk}
\alias{pamk}
\title{Clusters rare variants using fpc::pamk}
\usage{
  pamk(rv.dat, raw.dat)
}
\arguments{
  \item{rv.dat}{rare variant data frame}

  \item{raw.dat}{genotype/phenotype data frame}

  \item{vars}{list of strings identifying the cluster
  columns}

  \item{key}{string identifying the label column}
}
\value{
  cluster object containing\cr \tabular{ll}{ snp.clusters:
  \tab list containing SNPs grouped by cluster\cr
  cluster.info: \tab data frame containing information on
  each cluster\cr rv.dat: \tab rare variant data frame with
  cluster assignment added }
}
\description{
  Wrapper for fpc::pamk standardizing input and output\cr
  for rvclust.
}
\examples{
init(rv.dat,raw.dat)
init(rv.dat,raw.dat,vars=c("POS","STATE"))
init(rv.dat,raw.dat,key="SNP")
init(rv.dat,raw.dat,vars=c("POS","STATE"),key="SNP")
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{init}}

  \code{\link{annotate}}
}

