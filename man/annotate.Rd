\name{annotate}
\alias{annotate}
\title{Uses GWAR to annotate the genetic variants in an
rvclustobject.}
\usage{
  annotate(rv, annotations = NA, file = NA)
}
\arguments{
  \item{rv}{rvclustobject}
}
\value{
  annotated rvclustobject
}
\description{
  The Genome-Wide Annotation Repository (GWAR) is a\cr
  resource provided by the Bush Lab at Vanderbilt\cr
  University's that offers genome-wide annotations.\cr More
  information on GWAR can be found at:\cr
  \url{http://gwar.mc.vanderbilt.edu/}\cr \cr This method
  uses GWAR to annotate the variants\cr in an rvclustobject
  with the annotations requested in\cr the annotations
  list.\cr \cr Annotations will be appended to the ped
  object\cr
}
\author{
  R Michael Sivley \email{mike.sivley@vanderbilt.edu}
}
\seealso{
  \code{\link{rvclust}}
}

