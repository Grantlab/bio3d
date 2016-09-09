\name{is.mol2}
\alias{is.mol2}
\title{ Is an Object of Class \sQuote{mol2}?  }
\description{
  Checks whether its argument is an object of class \sQuote{mol2}.
}
\usage{
  is.mol2(x)
}
\arguments{
  \item{x}{ an R object. }
}
\details{
  Tests if the object \sQuote{x} is of class \sQuote{mol2}
  (\code{is.mol2}), i.e. if \sQuote{x} has a
  \dQuote{class} attribute equal to \code{mol2}.
}
\value{
  TRUE if x is an object of class \sQuote{mol2} and FALSE otherwise
}

\seealso{
  \code{\link{read.mol2}}
}
\examples{
# Read a PDB file
mol <- read.mol2( system.file("examples/aspirin.mol2", package="bio3d") )
is.mol2(mol)
}
\keyword{ classes }
