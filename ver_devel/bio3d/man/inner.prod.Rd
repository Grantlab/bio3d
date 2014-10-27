\name{inner.prod}
\alias{inner.prod}
\title{ Mass-weighted Inner Product }
\description{
  Inner product of vectors (mass-weighted if requested).
}
\usage{
inner.prod(x, y, mass=NULL)
}
\arguments{
  \item{x}{ a numeric vector or matrix. }
  \item{y}{ a numeric vector or matrix. }
  \item{mass}{ a numeric vector containing the atomic masses for
    weighting. }
}
\details{
  This function calculates the inner product between two vectors, or
  alternatively, the column-wise vector elements of matrices. If atomic
  masses are provided, the dot products will be mass-weighted.
  
  See examples for more details.
}
\value{
  Returns the inner product(s).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma} }, \code{\link{normalize.vector} } }
\examples{
## Matrix operations
x <- 1:3
y <- diag(x)
z <- matrix(1:9, ncol = 3, nrow = 3)

inner.prod(x,y)
inner.prod(y,z)


## Application to normal modes
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Check for orthogonality 
inner.prod(modes$U[,7], modes$U[,8])

}
\keyword{ utilities }
