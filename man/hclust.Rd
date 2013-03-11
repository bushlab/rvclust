\name{hclust}
\alias{hclust}
\title{Clusters rare variants using stats::hclust}
\usage{
  hclust(rv.dat, raw.dat, vars = c("POS", "STATE"),
    key = "SNP", objfun = power, ...)
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
  Wrapper for stats::hclust standardizing input and
  output\cr for rvclust. Frontier is computed using a
  statistical power proxy. User may specifiy an alternative
  objective function. If an objective function is
  specified, it should accept a vector of values as its
  first parameter. Any arguments to rvclust::hclust
  following the objective function will be passed as
  additional parameters to the function.
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{init}}

  \code{\link{annotate}}
}

