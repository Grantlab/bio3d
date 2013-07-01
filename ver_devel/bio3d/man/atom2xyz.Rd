\name{atom2xyz}
\alias{atom2xyz}
\alias{xyz2atom}
\title{ Convert Between Atom and xyz Indices }
\description{
  Basic functions to convert between xyz and their corresponding atom indices.
}
\usage{
atom2xyz(num)
xyz2atom(xyz.ind)
}
\arguments{
  \item{num}{ a numeric vector of atom indices. }
  \item{xyz.ind}{ a numeric vector of xyz indices. }
}
\value{
  A numeric vector of either xyz or atom indices.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{atom.select}}, \code{\link{read.pdb}} }
\examples{
xyz.ind <- atom2xyz(c(1,10,15))
xyz2atom( xyz.ind )
}
\keyword{ utilities }
