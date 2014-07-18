\name{mktrj.enma}
\alias{mktrj.enma}
\title{ Ensemble NMA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given normal mode
  vector.
}
\usage{
\method{mktrj}{enma}(enma = NULL, pdbs = NULL, ind = 1, mode = 1, mag = 10,
                     step = 1.25, file = NULL,  ...) 
}
\arguments{
  \item{enma}{ a list object of class \code{"enma"} (obtained with
    \code{\link{nma.pdbs}}).}
  \item{pdbs}{ a list object of class \code{"3dalign"} (obtained with
    \code{\link{pdbaln}} or \code{\link{read.fasta.pdb}}) which
    corresponds to the \code{"enma"} object.}
  \item{ind}{ structure id poiting to which structure in the \code{enma}
    object to make the trajectory of. }
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
\seealso{ \code{\link{nma}}, \code{\link{view.modes}}, \code{\link{mktrj.nma}},
  \code{\link{nma.pdbs}}, \code{\link{pdbaln}} }
\keyword{ utilities }
