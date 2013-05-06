\name{mktrj.pca}
\alias{mktrj.pca}
\title{ PCA Atomic Displacement Trajectory }
\description{
  Make a trajectory of atomic displacments along a given principal component.
}
\usage{
mktrj.pca(pca = NULL, pc = 1, mag = 1, step = 0.125, file = NULL,  ...)
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
\note{ The molecular graphics software VMD is useful for viewing
  trajectories see: \cr
  \url{http://www.ks.uiuc.edu/Research/vmd/}. }
\seealso{ \code{\link{pca.xyz}} }
\examples{
data(kinesin) # see pca.xyz for details
attach(kinesin, warn.conflicts=FALSE)

# Ignore gap containing positions
cut.seqs <- which(pdbs$id \%in\% c("d1n6mb_","d1ry6a_"))
gaps <- gap.inspect(pdbs$ali[-cut.seqs,])

# Write PC trajectory
a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb",
               resno = pdbs$resno[1, gaps$f.inds],
               resid = aa123(pdbs$ali[1, gaps$f.inds]) )

b <- mktrj.pca(pc.xray, pc=2, file="pc2.pdb",
               resno = pdbs$resno[1, gaps$f.inds],
               resid = aa123(pdbs$ali[1, gaps$f.inds]) )

c <- mktrj.pca(pc.xray, pc=3, file="pc3.pdb",
               resno = pdbs$resno[1, gaps$f.inds],
               resid = aa123(pdbs$ali[1, gaps$f.inds]) )

}
\keyword{ utilities }
