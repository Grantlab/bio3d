\name{inner.prod}
\alias{inner.prod}
\title{ Mass-weighted Inner Product }
\description{
  Inner product of vectors (mass-weighted if requested).
}
\usage{
inner.prod(a, b, mass=NULL)
}
\arguments{
  \item{a}{ a numeric vector or matrix. }
  \item{b}{ a numeric vector or matrix. }
  \item{mass}{ a numeric vector containing the atomic masses for
    weighting. }
}
\details{
  This function calculates the inner product between two vectors, or
  alternatively, the column-wise vector elements of matrices. If atomic
  masses are provided, the dot product will be mass-weighted.
  
  See examples for more details.
}
\value{
  Returns the normalized vector (or matrix).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma} }, \code{\link{normalize.vector} } }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Orthonormal mode vectors ?
ip <- inner.prod(normalize.vector(modes$U[,7]),
                 normalize.vector(modes$U[,8]))

ip <- inner.prod(normalize.vector(modes$U[,7:10]),
                 normalize.vector(modes$U[,7]))
}
\keyword{ utilities }
