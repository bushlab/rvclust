\name{burden.test}
\alias{burden.test}
\title{Collapses features into composite feature using the OR operator}
\usage{
  burden.test(rv)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  rvclustobject with collapsed data (composite features)
}
\description{
  Performs a standard burden test, which collapses the
  features in a bin by taking the SUM of all features. For
  each subject, the composite feature is represented by the
  burden of events occurring in the bin.
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}

  \code{\link{pamk}}
}

