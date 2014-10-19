\name{mktrj.nma}
\alias{mktrj.nma}
\title{ NMA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given normal mode
  vector.
}
\usage{
\method{mktrj}{nma}(nma = NULL, mode = 7, mag = 10, step = 1.25, file = NULL,  ...)
}
\arguments{
  \item{nma}{ a list object of class \code{"nma"} (obtained with
    \code{\link{nma}}).}
  \item{mode}{ the mode number along which displacements should be made.}
  \item{mag}{ a magnification factor for scaling the displacements. }
  \item{step}{ the step size by which to increment along the mode. }
  \item{file}{ a character vector giving the output PDB file name. }
  \item{...}{ extra arguments to be passed to the function write.pdb. }
}
\details{
  Trajectory frames are built from reconstructed Cartesian coordinates
  produced by interpolating from the structure along a given mode
  vector, in increments of \code{step}.

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
\note{ Molecular graphics software such as VMD or PyMOL is useful
  for viewing trajectories see e.g: \cr
  \url{http://www.ks.uiuc.edu/Research/vmd/}. }
\seealso{ \code{\link{nma}}, \code{\link{view.modes}}. }
\examples{
## Fetch stucture
pdb <- read.pdb("1hel")

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Visualize modes
outfile = file.path(tempdir(), "mode_7.pdb")
m7 <- mktrj(modes, mode=7, file = outfile)
outfile

}
\keyword{ utilities }
