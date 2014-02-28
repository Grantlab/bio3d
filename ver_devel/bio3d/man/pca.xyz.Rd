\name{pca.xyz}
\alias{pca.xyz}
\title{ Principal Component Analysis }
\description{
   Performs principal components analysis (PCA) on a \code{xyz}
   numeric data matrix.
 }
\usage{
pca.xyz(xyz, subset = rep(TRUE, nrow(as.matrix(xyz))), use.svd = FALSE)
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
}
\value{
  Returns a list with the following components:
  \item{L }{eigenvalues.}
  \item{U }{eigenvectors (i.e. the x, y, and z variable loadings).}
  \item{z }{scores of the supplied \code{xyz} on the pcs.}
  \item{au }{atom-wise loadings (i.e. xyz normalised eigenvectors).}
  \item{sdev }{the standard deviations of the pcs.}
  \item{mean }{the means that were subtracted.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{plot.pca}}, \code{\link{mktrj.pca}},
  \code{\link{pca.tor}}, \code{\link{pca.project}} }
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

#-- OR Read previously saved transducin data
data(transducin)
attach(transducin, warn.conflicts=FALSE)

# Previously fitted coordinates based on sub 1.0A^3 core
xyz <- pdbs$xyz

# Alternatively fit structures onto sub 0.5A^3 core
#xyz <- fit.xyz( fixed = pdbs$xyz[1,],
#                mobile = pdbs,
#                fixed.inds  = core$c0.5A.xyz,
#                mobile.inds = core$c0.5A.xyz) 
                
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

#-- Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])

# Plot results (conformer plots & scree plot)
plot(pc.xray, col=annotation[, "color"])

## Plot atom wise loadings
plot.bio3d(pc.xray$au[,1], ylab="PC1 (A)")

#plot.bio3d(gaps.res$f.ind, pc.xray$au[,1],
#           xlab="Alignment Position", ylab="PC1 (A)")

## Plot loadings in relation to reference structure "1TAG_A"
pdb <- read.pdb("1tag")

ind <- grep("1TAG", pdbs$id)
res.ref <- which(!is.gap(pdbs$ali[ind,]))
res.ind <- which(res.ref \%in\% gaps.res$f.ind)
op <- par(no.readonly=TRUE)
par(mfrow = c(3, 1), cex = 0.6, mar = c(3, 4, 1, 1))
plot.bio3d(res.ind, pc.xray$au[,1], sse=pdb, ylab="PC1 (A)")
plot.bio3d(res.ind, pc.xray$au[,2], sse=pdb, ylab="PC2 (A)")
plot.bio3d(res.ind, pc.xray$au[,3], sse=pdb, ylab="PC3 (A)")
par(op)

\dontrun{
# Write PC trajectory
a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )

b <- mktrj.pca(pc.xray, pc=2, file="pc2.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )

c <- mktrj.pca(pc.xray, pc=3, file="pc3.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )
}

detach(transducin)
}
\keyword{ utilities }
\keyword{ multivariate }

