\name{pairwise}
\alias{pairwise}
\title{ Pair Indices }
\description{
  A utility function to determine indices for pairwise comparisons.
}
\usage{
pairwise(N)
}
\arguments{
  \item{N}{ a single numeric value representing the total number of
    things to undergo pairwise comparison. }
}
\value{
  Returns a two column numeric matrix giving the indices for all
  pairs.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{seqidentity}} }
\examples{
  pairwise(3)
  pairwise(20)
}
\keyword{ utilities }
