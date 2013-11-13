\name{dccm.nma}
\alias{dccm.nma}
\title{ Dynamic Cross-Correlation from Normal Modes Analysis }
\description{
  Calculate the cross-correlation matrix from Normal Modes Analysis.
}
\usage{
\method{dccm}{nma}(x, nmodes = NULL, ncore = NULL, \dots)
}
\arguments{
  \item{x}{ an object of class \code{nma} as obtained from
    function \code{nma}. }
  \item{nmodes}{ numerical, number of modes to consider. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
  \item{\dots}{ additional arguments ? }
}
\details{
  This function calculates the cross-correlation matrix from Normal
  Modes Analysis (NMA) obtained from \code{nma} of a protein
  structure.  It returns a matrix of residue-wise cross-correlations
  whose elements, Cij, may be displayed in a graphical
  representation frequently termed a dynamical cross-correlation
  map, or DCCM.

  If Cij = 1 the fluctuations of residues i and j are completely
  correlated (same period and same phase), if Cij = -1 the
  fluctuations of residues i and j are completely anticorrelated (same
  period and opposite phase), and if Cij = 0 the fluctuations of i
  and j are not correlated.
}
\value{
  Returns a cross-correlation matrix.
}
\references{
  Wynsberghe. A.W.V, Cui, Q. \emph{Structure} \bold{14}, 1647--1653.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{plot.dccm}} }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate normal modes
modes <- nma(pdb)

## Calculate correlation matrix
cm <- dccm.nma(modes)

## Plot correlation map
plot(cm, sse = pdb, contour = F, col.regions = bwr.colors(20),
     at = seq(-1, 1, 0.1))
}
\keyword{ analysis }
