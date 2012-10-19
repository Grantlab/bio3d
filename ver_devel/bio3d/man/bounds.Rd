\name{bounds}
\alias{bounds}
\title{ Bounds of a Numeric Vector }
\description{
  Find the \sQuote{bounds} (i.e. start, end and length) of consecutive
  numbers within a larger set of numbers in a given vector.
}
\usage{
  bounds(nums, dup.inds=FALSE, pre.sort=TRUE)
}
\arguments{
  \item{nums}{ a numeric vector. }
  \item{dup.inds}{ logical, if TRUE the bounds of consecutive duplicated
    elements are returned. }
  \item{pre.sort}{ logical, if TRUE the input vector is ordered prior to
    bounds determination. }
}
\details{
   This is a simple utility function useful for summarizing the contents
   of a numeric vector. For example: find the start position, end
   position and lengths of secondary structure elements given a vector
   of residue numbers obtained from a DSSP secondary structure
   prediction.

   By setting \sQuote{dup.inds} to TRUE then the indices of the first
   (start) and last (end) duplicated elements of the vector are
   returned. For example: find the indices of atoms belonging to a particular
   residue given a vector of residue numbers (see below).
}
\value{
  Returns a three column matrix listing starts, ends and lengths.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22},
  2695--2696.
}
\author{ Barry Grant }
\examples{
 test <- c(seq(1,5,1),8,seq(10,15,1))
 bounds(test)

 test <- rep(c(1,2,4), times=c(2,3,4))
 bounds(test, dup.ind=TRUE)
}
\keyword{ utilities }
