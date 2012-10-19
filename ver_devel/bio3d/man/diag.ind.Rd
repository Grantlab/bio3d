\name{diag.ind}
\alias{diag.ind}
\title{ Diagonal Indices of a Matrix }
\description{
  Returns a matrix of logicals the same size of a given matrix with
  entries 'TRUE' in the upper triangle close to the diagonal.
}
\usage{
diag.ind(x, n = 1, diag = TRUE)
}
\arguments{
  \item{x}{ a matrix. }
  \item{n}{ the number of elements from the diagonal to include. }
  \item{diag}{ logical. Should the diagonal be included? }
}
\details{
  Basic function useful for masking elements close to the diagonal of a
  given matrix.
}
\value{
  Returns a matrix of logicals the same size of a given matrix with
  entries 'TRUE' in the upper triangle close to the diagonal.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{diag}},  \code{\link{lower.tri}},
  \code{\link{upper.tri}}, \code{\link{matrix}} }
\examples{
diag.ind( matrix(,ncol=5,nrow=5), n=3 )
}
\keyword{ utilities }
