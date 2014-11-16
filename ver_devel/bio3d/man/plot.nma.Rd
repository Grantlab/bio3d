\name{plot.nma}
\alias{plot.nma}
\title{ Plot NMA Results }
\description{
  Produces eigenvalue/frequency spectrum plots and an atomic fluctuations plot.
}
\usage{
\method{plot}{nma}(x, pch = 16, col = par("col"), cex=0.8, mar=c(6, 4, 2, 2),...)
}
\arguments{
  \item{x}{ the results of normal modes analysis obtained with
    \code{\link{nma}}. }
  \item{pch}{ a vector of plotting characters or symbols: see \sQuote{points}. }
  \item{col}{ a character vector of plotting colors.  }
  \item{cex}{ a numerical single element vector giving the amount by
    which plotting text and symbols should be magnified relative to the
    default. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{\dots}{ extra plotting arguments passed to \code{plot.bio3d}
    that effect the atomic fluctuations plot only. }
}
\details{
  \code{plot.nma} produces an eigenvalue (or frequency) spectrum plot
  together with a plot of the atomic fluctuations. 
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{plot.bio3d}} }
\examples{
## Fetch structure
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate modes
modes <- nma(pdb)

plot(modes, sse=pdb)

}
\keyword{ hplot }
