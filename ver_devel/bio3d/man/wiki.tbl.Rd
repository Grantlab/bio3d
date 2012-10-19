\name{wiki.tbl}
\alias{wiki.tbl}
\title{ Wiki Table }
\description{
  Print a matrix as a Wiki formatted table.
}
\usage{
wiki.tbl(mat)
}
\arguments{
  \item{mat}{ a matrix for printing. }
}
\details{
  This function prints a basic wiki formatted table from a given matrix \sQuote{mat}.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\examples{
data(USArrests)
wiki.tbl(USArrests)
# paste output into your PBWIKI
}
\keyword{ utilities }
