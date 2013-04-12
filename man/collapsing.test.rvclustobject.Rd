\name{collapsing.test.rvclustobject}
\alias{collapsing.test.rvclustobject}
\title{Collapses features into composite feature using the OR operator}
\usage{
  \method{collapsing.test}{rvclustobject} (rv)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  rvclustobject with collapsed data (composite features)
}
\description{
  Performs a standard collapsing test, which collapses the
  features in a bin by taking the OR of all features. For
  each subject, the composite feature is represented as
  positive if any of the individual events have occurred.
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclustobject}}

  \code{\link{annotate}}

  \code{\link{pamk}}
}

