\name{mktrj}
\alias{mktrj}
\title{ PCA / NMA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given principal
  component / normal mode.
}
\usage{
mktrj(...)
}
\arguments{
  \item{\dots}{ arguments passed to the methods
    \code{mktrj.pca}, \code{mktrj.nma}, or \code{mktrj.enma}. }
}
\details{
  \code{mktrj} is a generic function calling the corresponding function
  determined by the class of the input argument \code{x}. Use
  \code{methods("mktrj")} to get all the methods for \code{mktrj}
  generic:
  
  \code{\link{mktrj.pca}} will be used when \code{x} is an object of 
  \code{"pca"}.
  
  \code{\link{mktrj.nma}} will be used when \code{x} is an object of 
  \code{"nma"}.
  
  See examples for each corresponding function for more details.
}
\note{ Molecular graphics software such as VMD or PyMOL is useful
  for viewing trajectories see e.g: \cr
  \url{http://www.ks.uiuc.edu/Research/vmd/}. }
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{ \code{\link{mktrj.pca}}, \code{\link{mktrj.nma}},
  \code{\link{pca.xyz}}, \code{\link{nma}}, \code{\link{view.modes}}. }
\examples{
data(transducin)
attach(transducin, warn.conflicts=FALSE)

# Ignore gap containing positions
gaps.pos <- gap.inspect(pdbs$xyz)

# PCA
pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])

# Write PC trajectory of pc=1
a <- mktrj(pc.xray)

detach(transducin)
}
\keyword{ utilities }
