\name{init}
\alias{init}
\title{Initializes data and covariates for rare variant analysis}
\usage{
  init(data.path, data.fname, cov.path = NA,
    cov.fname = NA)
}
\arguments{
  \item{data.path}{data directory}

  \item{data.fname}{data file name}

  \item{cov.path}{covariate file directory}

  \item{cov.fname}{covariate file name [optional]}
}
\value{
  rvclust object containing\cr \tabular{ll}{ rv.dat: \tab
  rare variant data frame\cr raw.dat: \tab
  genotype/phenotype data frame\cr cov.dat: \tab covariate
  data frame [or NA] }
}
\description{
  Utility to load data, phenotype, and an optional
  covariates\cr file into data frames formatted for rare
  variant clustering\cr and analysis with rvclust.
}
\note{
  Data files should be in PLINK ped/map format

  Covariates file should have one column for subject

  Do not include ped/map extension

  Do include covariate extension IDs followed by one column
  per covariate
}
\examples{
init('/path/to','data')
init('/path/to','data','/path/to','covariates.csv')
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}

