\name{hpower}
\alias{hpower}
\title{Clusters rare variants using rvclust custom hierarchical method}
\usage{
  hpower(rv.dat, raw.dat, vars = c("POS", "STATE"),
    key = "SNP")
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
  cluster assignment added\cr }
}
\description{
  Top down hierarchical cluster where each level is
  partitioned\cr into two clusters using fpc::pamk across
  the specified columns.\cr The fitness function, which
  decides whether to expand a cluster\cr is a statistical
  power proxy provided by Dr. Eli Stahl, Emory University.
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

