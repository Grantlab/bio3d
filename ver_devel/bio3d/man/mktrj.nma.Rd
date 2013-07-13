\name{mktrj.nma}
\alias{mktrj.nma}
\title{ NMA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given normal mode
  vector.
}
\usage{
mktrj.nma(nma = NULL, mode = 7, mag = 100, step = 12.5, file = NULL,  ...)
}
\arguments{
  \item{nma}{ a list object of class \code{"nma"} (obtained with
    \code{\link{nma}}).}
  \item{mode}{ the mode number along which displacements should be made.}
  \item{mag}{ a magnification factor for scaling the displacements. }
  \item{step}{ the step size by which to increment along the \code{pc}. }
  \item{file}{ a character vector giving the output PDB file name. }
  \item{...}{ extra arguments to be passed to the function write.pdb. }
}
\details{
  Trajectory frames are built from reconstructed Cartesian coordinates
  produced by interpolating from the mean structure along a given
  \code{normal mode}, in increments of \code{step}.

  An optional magnification factor can be used to amplify
  displacements.  This involves scaling the mode vector by
  \code{mag}-times.
}
\value{
  Returns a numeric matrix of interpolated coordinates with a row per structure.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\note{ The molecular graphics software VMD is useful for viewing
  trajectories see: \cr
  \url{http://www.ks.uiuc.edu/Research/vmd/}. }
\seealso{ \code{\link{nma} } }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Visualize modes
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")
}
\keyword{ utilities }
