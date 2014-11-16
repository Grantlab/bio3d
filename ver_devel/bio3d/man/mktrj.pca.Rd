\name{mktrj.pca}
\alias{mktrj.pca}
\title{ PCA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given principal component.
}
\usage{
\method{mktrj}{pca}(pca = NULL, pc = 1, mag = 1, step = 0.125, file = NULL,  ...)
}
\arguments{
  \item{pca}{ a list object of class \code{"pca"} (obtained with
    \code{\link{pca.xyz}}).}
  \item{pc}{ the PC number along which displacements should be made.}
  \item{mag}{ a magnification factor for scaling the displacements. }
  \item{step}{ the step size by which to increment along the \code{pc}. }
  \item{file}{ a character vector giving the output PDB file name. }
  \item{...}{ extra arguments to be passed to the function write.pdb. }  
}
\details{
  Trajectory frames are built from reconstructed Cartesian coordinates
  produced by interpolating from the mean structure along a given
  \code{pc}, in increments of \code{step}.

  An optional magnification factor can be used to amplify
  displacements.  This involves scaling by \code{mag}-times the standard
  deviation of the conformer distribution along the given \code{pc}
  (i.e. the square root of the associated eigenvalue).
}
\value{
  Returns a numeric matrix of interpolated coordinates with a row per structure.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ Molecular graphics software such as VMD or PyMOL is useful
  for viewing trajectories see e.g: \cr
  \url{http://www.ks.uiuc.edu/Research/vmd/}. }
\seealso{ \code{\link{pca.xyz}}, \code{\link{view.modes}}. }
\examples{
\donttest{
data(transducin)
attach(transducin, warn.conflicts=FALSE)

# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

# PCA
pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])

# Write PC trajectory
path = tempdir()
a <- mktrj(pc.xray, pc=1, file=file.path(path, "pc1.pdb"),
           resno = pdbs$resno[1, gaps.res$f.inds],
           resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )

b <- mktrj(pc.xray, pc=2, file=file.path(path, "pc2.pdb"),
           resno = pdbs$resno[1, gaps.res$f.inds],
           resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )

c <- mktrj(pc.xray, pc=3, file=file.path(path, "pc3.pdb"),
           resno = pdbs$resno[1, gaps.res$f.inds],
           resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )

# Output files
file.path(path, paste("pc", 1:3, ".pdb", sep=""))

detach(transducin)
}
}
\keyword{ utilities }
