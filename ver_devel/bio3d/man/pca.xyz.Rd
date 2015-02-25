\name{pca.xyz}
\alias{pca.xyz}
\alias{print.pca}
\title{ Principal Component Analysis }
\description{
   Performs principal components analysis (PCA) on a \code{xyz}
   numeric data matrix.
 }
\usage{
\method{pca}{xyz}(xyz, subset = rep(TRUE, nrow(as.matrix(xyz))),
                 use.svd = FALSE, rm.gaps=FALSE, mass = NULL, \dots)

\method{print}{pca}(x, nmodes=6, \dots)
}
\arguments{
  \item{xyz}{ numeric matrix of Cartesian coordinates with a row per
    structure. }
  \item{subset}{ an optional vector of numeric indices that selects a
    subset of rows (e.g. experimental structures vs molecular dynamics
    trajectory structures) from the full \code{xyz} matrix. Note: the
    full \code{xyz} is projected onto this subspace.}
  \item{use.svd}{ logical, if TRUE singular value decomposition (SVD) is 
    called instead of eigenvalue decomposition. }
  \item{rm.gaps}{ logical, if TRUE gap positions (with missing
    coordinate data in any input structure) are removed before
    calculation. This is equivalent to removing NA cols from xyz. }
  \item{x}{ an object of class \code{pca}, as obtained from function
    \code{pca.xyz}. }
  \item{nmodes}{ numeric, number of modes to be printed. }
  \item{mass}{ a \sQuote{pdb} object or numeric vector of residue/atom masses. 
     By default (\code{mass=NULL}), mass is ignored. If provided 
     with a \sQuote{pdb} object, masses of all amino acids obtained from 
     \code{\link{aa2mass}} are used. }
  \item{\dots}{ additional arguments to \code{\link{fit.xyz}} 
     (for \code{pca.xyz}) or to \code{print} (for \code{print.pca}). }
}
\note{
     If \code{mass} is provided, mass weighted coordinates will be considered, 
     and iteration of fitting onto the mean structure is performed internally. 
     The extra fitting process is to remove external translation and rotation 
     of the whole system. With this option, a direct comparison can be made 
     between PCs from \code{\link{pca.xyz}} and vibrational modes from 
     \code{\link{nma.pdb}}, with the fact that 
                           \deqn{A=k_BTF^{-1}}{A=k[B]TF^-1}, 
     where \eqn{A} is the variance-covariance matrix, \eqn{F} the Hessian 
     matrix, \eqn{k_B}{k[B]} the Boltzmann's constant, and \eqn{T} the 
     temperature.
}
\value{
  Returns a list with the following components:
  \item{L }{eigenvalues.}
  \item{U }{eigenvectors (i.e. the x, y, and z variable loadings).}
  \item{z }{scores of the supplied \code{xyz} on the pcs.}
  \item{au }{atom-wise loadings (i.e. xyz normalized eigenvectors).}
  \item{sdev }{the standard deviations of the pcs.}
  \item{mean }{the means that were subtracted.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{
  \code{\link{pca}}, \code{\link{pca.pdbs}},
  \code{\link{plot.pca}}, \code{\link{mktrj.pca}},
  \code{\link{pca.tor}}, \code{\link{project.pca}} }
\examples{

\dontrun{
#-- Read transducin alignment and structures
aln <- read.fasta(system.file("examples/transducin.fa",package="bio3d"))
pdbs <- read.fasta.pdb(aln)

# Find core
core <- core.find(pdbs, 
                  #write.pdbs = TRUE,
                  verbose=TRUE)

rm(list=c("pdbs", "core"))
}

#-- OR for demo purposes just read previously saved transducin data
attach(transducin)

# Previously fitted coordinates based on sub 1.0A^3 core. See core.find() function.
xyz <- pdbs$xyz
                
#-- Do PCA ignoring gap containing positions
pc.xray <- pca.xyz(xyz, rm.gaps=TRUE)

# Plot results (conformer plots & scree plot overview)
plot(pc.xray, col=annotation[, "color"])

# Plot a single conformer plot of PC1 v PC2
plot(pc.xray, pc.axes=1:2, col=annotation[, "color"])

## Plot atom wise loadings
plot.bio3d(pc.xray$au[,1], ylab="PC1 (A)")

\donttest{
## Plot loadings in relation to reference structure 1TAG
pdb <- read.pdb("1tag")
ind <- grep("1TAG", pdbs$id)             ## location in alignment

resno <- pdbs$resno[ind, !is.gap(pdbs)]  ## non-gap residues
tpdb <- trim.pdb(pdb, resno=resno)

op <- par(no.readonly=TRUE)
par(mfrow = c(3, 1), cex = 0.6, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno, ylab="PC1 (A)", sse=tpdb)
plot.bio3d(pc.xray$au[,2], resno, ylab="PC2 (A)", sse=tpdb)
plot.bio3d(pc.xray$au[,3], resno, ylab="PC3 (A)", sse=tpdb)
par(op)
}

\dontrun{
# Write PC trajectory
resno = pdbs$resno[1, !is.gap(pdbs)]
resid = aa123(pdbs$ali[1, !is.gap(pdbs)])

a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb",
               resno=resno, resid=resid )

b <- mktrj.pca(pc.xray, pc=2, file="pc2.pdb",
               resno=resno, resid=resid )

c <- mktrj.pca(pc.xray, pc=3, file="pc3.pdb",
               resno=resno, resid=resid )
}

detach(transducin)
}
\keyword{ utilities }
\keyword{ multivariate }

