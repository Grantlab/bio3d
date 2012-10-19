\name{atom2xyz}
\alias{atom2xyz}
\title{ Convert Between Atom and xyz Indices }
\description{
  Basic function to return the Calpha xyz indices given their
  corresponding Calpha only atom indices.
}
\usage{
atom2xyz(num)
}
\arguments{
  \item{num}{ a numeric vector of atom indices. }
}
\value{
  A numeric vector of xyz indices.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{atom.select}}, \code{\link{read.pdb}} }
\examples{
xyz.ind <- atom2xyz(c(1,10,15))

## For xyz2atom:
unique( ceiling(xyz.ind/3) )
}
\keyword{ utilities }
