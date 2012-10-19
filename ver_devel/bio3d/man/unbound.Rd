\name{unbound}
\alias{unbound}
\title{ Sequence Generation from a Bounds Vector }
\description{
  Generate a sequence of consecutive numbers from a \code{\link{bounds}} vector.
}
\usage{
unbound(start, end)
}
\arguments{
  \item{start}{ vector of starting values, such as that obtained from \code{\link{bounds}}. }
  \item{end}{ vector of (maximal) end values, such as that obtained from \code{\link{bounds}}. }
}
\details{
   This is a simple utility function that does the opposite of the
   \code{\link{bounds}} function.
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
unbound(b[,"start"],b[,"end"])

}
\keyword{ utilities }
