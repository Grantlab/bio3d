\name{is.pdb}
\alias{is.pdb}
\alias{is.pdbs}
\title{ Is an Object of Class \sQuote{pdb(s)}?  }
\description{
  Checks whether its argument is an object of class \sQuote{pdb} or \sQuote{pdbs}.
}
\usage{
  is.pdb(x)
  is.pdbs(x)
}
\arguments{
  \item{x}{ an R object. }
}
\details{
  Tests if the object \sQuote{x} is of class \sQuote{pdb}
  (\code{is.pdb}) or \sQuote{pdbs} (\code{is.pdbs}), i.e. if \sQuote{x} has a
  \dQuote{class} attribute equal to \code{pdb} or \code{pdbs}.
}
\value{
  TRUE if x is an object of class \sQuote{pdb(s)} and FALSE otherwise
}

\seealso{
  \code{\link{read.pdb}}, \code{\link{read.fasta.pdb}},
  \code{\link{pdbaln}}
}
\examples{
# Read a PDB file
pdb <- read.pdb("1BG2")
is.pdb(pdb)
}
\keyword{ classes }
