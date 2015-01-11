\name{rle2}
\title{Run Length Encoding with Indices}
\alias{rle2}
\alias{print.rle2}
\concept{runs}
\description{
  Compute the lengths, values and indices of runs of equal values in a
  vector. This is a modifed version of base function \code{rle()}.
}
\usage{
rle2(x)

\method{print}{rle2}(x, digits = getOption("digits"), prefix = "", \dots)
}
\arguments{
  \item{x}{an atomic vector for \code{rle()}; an object of class
    \code{"rle"} for \code{inverse.rle()}.}
  \item{\dots}{further arguments; ignored here.}
  \item{digits}{number of significant digits for printing, see
    \code{\link{print.default}}.}
  \item{prefix}{character string, prepended to each printed line.}
}
\details{
  Missing values are regarded as unequal to the previous value, even if
  that is also missing.

  \code{inverse.rle()} is the inverse function of \code{rle2()} and \code{rle()},
  reconstructing \code{x} from the runs.
}
\value{
  \code{rle()} returns an object of class \code{"rle"} which is a list
  with components:
  \item{lengths}{an integer vector containing the length of each run.}
  \item{values}{a vector of the same length as \code{lengths} with the
    corresponding values.}
}
\examples{
x <- rev(rep(6:10, 1:5))
rle(x)
## lengths [1:5]  5 4 3 2 1
## values  [1:5] 10 9 8 7 6
rle2(x)
## lengths: int [1:5] 5 4 3 2 1
## values : int [1:5] 10 9 8 7 6
## indices: int [1:5] 5 9 12 14 15
}
\keyword{manip}
