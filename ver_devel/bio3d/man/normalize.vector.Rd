\name{normalize.vector}
\alias{normalize.vector}
\title{ Mass-Weighted Normalized Vector }
\description{
  Normalizes a vector (mass-weighted if requested).
}
\usage{
normalize.vector(v, mass=NULL)
}
\arguments{
  \item{v}{ a numeric vector or matrix to be normalized. }
  \item{mass}{ a numeric vector containing the atomic masses for
    weighting. }
}
\details{
  This function normalizes a vector, or alternatively, the column-wise
  vector elements of a matrix. If atomic masses are provided the vector
  is mass-weigthed. 
  
  See examples for more details.
}
\value{
  Returns the normalized vector (or matrix).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma} }, \code{\link{inner.prod} } }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Returns a vector
nv <- normalize.vector(modes$modes[,7], mass=modes$mass)

## Returns a matrix
nv <- normalize.vector(modes$modes[,7:10], mass=modes$mass)
}
\keyword{ utilities }
