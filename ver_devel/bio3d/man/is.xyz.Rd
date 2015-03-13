\name{is.xyz}
\alias{is.xyz}
\alias{as.xyz}
\title{ Is an Object of Class \sQuote{xyz}?  }
\description{
  Checks whether its argument is an object of class \sQuote{xyz}.
}
\usage{
  is.xyz(x)
  as.xyz(x)
}
\arguments{
  \item{x}{
    an R object to be tested
  }
}
\details{
  Tests if x is an object of class \sQuote{xyz}, i.e. if x has a
  \dQuote{class} attribute equal to \code{xyz}. 
}
\value{
  TRUE if x is an object of class \sQuote{xyz} and FALSE otherwise
}

\seealso{
  \code{\link{read.pdb}}, \code{\link{read.ncdf}},
  \code{\link{read.dcd}}, \code{\link{fit.xyz}}
}
\examples{
# Read a PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
is.xyz(pdb$xyz)

}
\keyword{ classes }
