\name{dist.xyz}
\alias{dist.xyz}
\title{  Calculate the Distances Between the Rows of Two Matrices }
\description{
  Compute the pairwise euclidean distances between the rows of two matrices.
}
\usage{
dist.xyz(a, b = NULL, all.pairs=TRUE, ncore=1, nseg.scale=1)
}
\arguments{
  \item{a}{ a \sQuote{xyz} object, numeric data matrix, or vector. }
  \item{b}{ an optional second \sQuote{xyz} object, data matrix, or vector. }
  \item{all.pairs}{  logical, if TRUE all pairwise distances between the
    rows of \sQuote{a} and all rows of \sQuote{b} are computed, if FALSE
    only the distances between coresponding rows of \sQuote{a} and
    \sQuote{b} are computed. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. } 
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
 This function returns a matrix of euclidean distances between each row
 of \sQuote{a} and all rows of \sQuote{b}.  Input vectors are
 coerced to three dimensional matrices (representing the Cartesian
 coordinates x, y and z) prior to distance computation.  If \sQuote{b} is
 not provided then the pairwise distances between all rows of \sQuote{a}
 are computed.
}
\value{
  Returns a matrix of pairwise euclidean distances between each row of
  \sQuote{a} and all rows of \sQuote{b}.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  This function will choke if \sQuote{b} has too many rows.
}
\seealso{ \code{\link{dm}}, \code{\link{dist}} }
\examples{
 dist.xyz( c(1,1,1, 3,3,3), c(3,3,3, 2,2,2, 1,1,1))
 dist.xyz( c(1,1,1, 3,3,3), c(3,3,3, 2,2,2, 1,1,1), all.pairs=FALSE)

}
\keyword{ utilities }
