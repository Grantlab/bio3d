\name{unbound}
\alias{unbound}
\title{ Sequence Generation from a Bounds Vector }
\description{
  Generate a sequence of consecutive numbers from a \code{\link{bounds}} vector.
}
\usage{
unbound(start, end = NULL)
}
\arguments{
  \item{start}{ vector of starting values, or a matrix containing starting and 
      end values such as that obtained from \code{\link{bounds}}. }
  \item{end}{ vector of (maximal) end values, such as that obtained from 
      \code{\link{bounds}}. }
}
\details{
   This is a simple utility function that does the opposite of the
   \code{\link{bounds}} function. If \code{start} is a vector, \code{end} must 
   be a vector having the same length as \code{start}. If \code{start} is a 
   matrix with column names contain 'start' and 'end', such as that returned 
   from \code{\link{bounds}}, \code{end} can be skipped and both starting and 
   end values will be extracted from \code{start}.
}
\value{
  Returns a numeric sequence vector.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\seealso{ \code{\link{bounds}} }
\author{ Barry Grant }
\examples{

test <- c(seq(1,5,1),8,seq(10,15,1))
b <- bounds(test)
unbound(b)

}
\keyword{ utilities }
