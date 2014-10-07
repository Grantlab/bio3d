\name{setup.ncore}
\alias{setup.ncore}
\title{ Setup for Running Bio3D Functions using Multiple CPU Cores }
\description{
  Internally used in parallelized Bio3D functions.
}
\usage{
  setup.ncore(ncore, bigmem = FALSE)
}
\arguments{
  \item{ncore}{
    User set (or default) value of \sQuote{ncore}.
  }
  \item{bigmem}{
    logical, if TRUE also check the availability of \sQuote{bigmemory} package.
  }
}
\details{
  Check packages and set correct value of \sQuote{ncore}.
}
\value{
  The actual value of \sQuote{ncore}.
}
\examples{
   setup.ncore(NULL)
   setup.ncore(1)
#   setup.ncore(2)
}
\keyword{ utilities }
