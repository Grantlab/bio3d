\name{store.atom}
\alias{store.atom}
\title{ Store all-atom data from a PDB object }
\description{
  Not intended for public usage
}
\usage{
store.atom(pdb)
}
\arguments{
  \item{pdb}{ A pdb object as obtained from read.pdb }
}
\details{
  This function was requested by a user and has not been extensively
  tested.  Hence it is not yet recommended for public usage.
}
\value{
  Returns a matrix of all-atom data
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ This function is still in development and is NOT part of the
  offical bio3d package
}
\seealso{ \code{\link{read.fasta.pdb}} }
\examples{
\dontrun{
pdb <- read.pdb( get.pdb("5p21", URLonly=TRUE) )
a <- store.atom(pdb)
a[,,1:2]
}
}
\keyword{ utilities }
