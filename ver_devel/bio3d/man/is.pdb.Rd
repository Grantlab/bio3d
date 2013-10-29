\name{is.pdb}
\alias{is.pdb}
\title{
  Is an Object of Class \sQuote{pdb}?
}
\description{
  Checks whether its argument is an object of class \sQuote{pdb}.
}
\usage{
  is.pdb(x)
}
\arguments{
  \item{x}{
    an R object to be tested
  }
}
\details{
  Tests if x is an object of class \sQuote{pdb}, i.e. if x has a \dQuote{class} attribute equal to \code{pdb}.
}
\value{
  TRUE if x is an object of class \sQuote{pdb} and FALSE otherwise
}

\seealso{
  \code{\link{read.pdb}}
}
\examples{
# Read a PDB file
pdb <- read.pdb("1BG2")
is.pdb(pdb)

}
\keyword{ classes }
