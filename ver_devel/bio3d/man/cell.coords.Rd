\encoding{UTF-8}
\name{cell.coords}
\alias{cell.coords}
\alias{cell.coords.default}
\title{Unit Cell Vectors Coordinates}
\description{
Computes the Cartesian coordinates of lattice vectors.
}
\usage{
cell.coords(\dots)

\method{cell.coords}{default}(x, digits = 3, \dots)
}
\arguments{
  \item{x}{a length 6 numeric vector containing the a, b, c, alpha, beta and gamma lattice parameters.}
  \item{digits}{an integer used to round the returned matrix}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
These functions compute a 3x3 matrix whose columns contrain the Cartesian coordinates of lattice vectors.
The 'a' and 'b' vectors are assumed to be respectively along the x-axis and in the xy-plane.
}
\value{
A 3x3 matrix containing the Cartesian coordinates of lattice vectors arranged by columns.
}

\author{
Julien Idé
}

\seealso{
To read unit cell parameters from trajectory files: \code{\link{read.ncdf}}, \code{\link{read.dcd}}
}
\examples{
cell <- c(a=20, b=20, c=20, alpha=90, beta=90, gamma=60)
cell.coords(cell)
cell.coords(cell, digits=6)
}
\keyword{ utilities }
