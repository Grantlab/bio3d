\name{mktrj}
\alias{mktrj}
\alias{mktrj.pca}
\alias{mktrj.nma}
\alias{mktrj.enma}
\title{ PCA / NMA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given principal
  component / normal mode.
}
\usage{
mktrj(...)

\method{mktrj}{pca}(pca = NULL, pc = 1, mag = 1, step = 0.125, file = NULL,  ...)

\method{mktrj}{nma}(nma = NULL, mode = 7, mag = 10, step = 1.25, file = NULL,  ...)

\method{mktrj}{enma}(enma = NULL, pdbs = NULL, s.inds = NULL, m.inds = NULL,
      mag = 10, step = 1.25, file = NULL, rock = TRUE, ncore = NULL, ...) 
}
\arguments{
  \item{pca}{ an object of class \code{"pca"} as obtained with
    function \code{\link{pca.xyz}} or \code{\link{pca}}.}
  \item{nma}{ an object of class \code{"nma"} as obtained with
    function \code{\link{nma.pdb}}.}
  \item{enma}{ an object of class \code{"enma"} as obtained with
    function \code{\link{nma.pdbs}}.}
  \item{pc}{ the PC number along which displacements should be made.}
  \item{mag}{ a magnification factor for scaling the displacements. }
  \item{step}{ the step size by which to increment along the pc/mode. }
  \item{file}{ a character vector giving the output PDB file name. }
  \item{mode}{ the mode number along which displacements should be made.}
  \item{pdbs}{ a list object of class \code{"pdbs"} (obtained with
    \code{\link{pdbaln}} or \code{\link{read.fasta.pdb}}) which
    corresponds to the \code{"enma"} object.}
  \item{s.inds}{ index or indices pointing to the structure(s) in the
    \code{enma} object for which the trajectory shall be generated. }
  \item{m.inds}{ the mode number(s) along which displacements should be
    made. }
  \item{rock}{ logical, if TRUE the trajectory rocks. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. } 
  \item{\dots}{ additional arguments passed to and from functions
    (e.g. to function \code{\link{write.pdb}}). }
}
\details{
  Trajectory frames are built from reconstructed Cartesian coordinates
  produced by interpolating from the mean structure along a given
  \code{pc} or \code{mode}, in increments of \code{step}.

  An optional magnification factor can be used to amplify
  displacements.  This involves scaling by \code{mag}-times the standard
  deviation of the conformer distribution along the given \code{pc}
  (i.e. the square root of the associated eigenvalue).
}
\note{ Molecular graphics software such as VMD or PyMOL is useful
  for viewing trajectories see e.g: \cr
  \url{http://www.ks.uiuc.edu/Research/vmd/}. }
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{
  \code{\link{pca}}, \code{\link{nma}},
  \code{\link{nma.pdbs}},
  \code{\link{view.modes}}.
}
\examples{
##- PCA example
data(transducin)
attach(transducin, warn.conflicts=FALSE)

# Calculate principal components
pc.xray <- pca(pdbs, fit=TRUE)

# Write PC trajectory of pc=1
outfile = tempfile()
a <- mktrj(pc.xray, file = outfile)
outfile
detach(transducin)


##- NMA example
## Fetch stucture
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Visualize modes
outfile = file.path(tempdir(), "mode_7.pdb")
mktrj(modes, mode=7, file = outfile)
outfile

}
\keyword{ utilities }
