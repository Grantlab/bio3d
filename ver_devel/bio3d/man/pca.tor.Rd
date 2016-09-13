\name{pca.tor}
\alias{pca.tor}
\title{ Principal Component Analysis }
\description{
   Performs principal components analysis (PCA) on torsion angle \code{data}.
 }
\usage{
\method{pca}{tor}(data, \dots)
}
\arguments{
  \item{data}{ numeric matrix of torsion angles with a row per structure. }
  \item{\dots}{ additional arguments passed to the method \code{pca.xyz}. }
}
\value{
  Returns a list with the following components:
  \item{L }{eigenvalues.}
  \item{U }{eigenvectors (i.e. the variable loadings).}
  \item{z.u }{scores of the supplied \code{data} on the pcs.}
  \item{sdev }{the standard deviations of the pcs.}
  \item{mean }{the means that were subtracted.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{  Barry Grant and Karim ElSawy }
\seealso{ \code{\link{torsion.xyz}}, \code{\link{plot.pca}},
  \code{\link{plot.pca.loadings}}, \code{\link{pca.xyz}} }
\examples{
##-- PCA on torsion data for multiple PDBs 
attach(kinesin)

gaps.pos <- gap.inspect(pdbs$xyz)
tor <- t(apply( pdbs$xyz[, gaps.pos$f.inds], 1, torsion.xyz, atm.inc=1))
pc.tor <- pca.tor(tor[,-c(1,218,219,220)])
#plot(pc.tor)
plot.pca.loadings(pc.tor)

detach(kinesin)

\dontrun{
##-- PCA on torsion data from an MD trajectory
trj <- read.dcd( system.file("examples/hivp.dcd", package="bio3d") )
tor <- t(apply(trj, 1, torsion.xyz, atm.inc=1))
gaps <- gap.inspect(tor)
pc.tor <- pca.tor(tor[,gaps$f.inds])
plot.pca.loadings(pc.tor)
}
}
\keyword{ utilities }
\keyword{ multivariate }

