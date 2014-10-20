\name{plot.rmsip}
\alias{plot.rmsip}
\title{ Plot RMSIP Results }
\description{
  Produces a heat plot of RMSIP (Root mean square inner product) for the
  visualization of modes similarity.
}
\usage{
\method{plot}{rmsip}(x, xlab = NULL, ylab = NULL, col = gray(50:0/50),
      zlim=c(0,1), \dots)
}
\arguments{
  \item{x}{ an object of class \code{rmsip}. }
  \item{xlab}{ a label for the x axis, defaults to \sQuote{a}. }
  \item{ylab}{ a label for the y axis, defaults to \sQuote{b}. }
  \item{col}{ a vector of colors for the RMSIP map (or overlap values). }
  \item{zlim}{ the minimum and maximum \sQuote{z} values for which
    colors should be plotted. }
  \item{\dots}{ additional arguments to function \code{image}. }
}
\details{
  \code{plot.rmsip} produces a color image with the function
  \code{image}. 
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{rmsip}}, \code{\link{overlap}},
  \code{\link{nma}}, \code{\link{image}}.  }
\examples{
## Read PDB structure
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Perform NMA
modes.a <- nma(pdb, ff="calpha")
modes.b <- nma(pdb, ff="anm")

## Calculate and plot RMSIP
r <- rmsip(modes.a, modes.b)
plot(r)
}
\keyword{ hplot }
