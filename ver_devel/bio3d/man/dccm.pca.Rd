\name{dccm.pca}
\alias{dccm.pca}
\title{ Dynamical Cross-Correlation Matrix from Principal Component Analysis }
\description{
  Calculate the cross-correlation matrix from principal component analysis (PCA).
}
\usage{
\method{dccm}{pca}(x, pc = NULL, method = c("pearson", "lmi"), ncore = NULL, \dots)
}
\arguments{
  \item{x}{ an object of class \code{pca} as obtained from
    function \code{pca.xyz}. }
  \item{pc}{ numerical, indices of PCs to be included in the calculation.
    If all negative, PCs complementary to \code{abs(pc)} are included. }
  \item{method}{ method to calculate the cross-correlation. Currently supports
    Pearson and linear mutual information (LMI). }
  \item{ncore}{ number of CPU cores used to do the calculation.
    By default (\code{ncore = NULL}), use all available cores detected. }
  \item{\dots}{ Additional arguments to be passed (currently ignored). }
}
\details{
  This function calculates the cross-correlation matrix from principal 
  component analysis (PCA) obtained from \code{pca.xyz} of a set of protein
  structures. It is an alternative way to calculate correlation in addition
  to the conventional way from xyz coordinates directly. But, in this new 
  way one can freely chooses the PCs to be included in the 
  calculation (e.g. for filtering out PCs with small eigenvalues).
}
\value{
  Returns a cross-correlation matrix with values in a range from -1 to 1 (Pearson)
  or from 0 to 1 (LMI).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao}
\seealso{ \code{\link{pca.xyz}}, \code{\link{plot.dccm}}, \code{\link{dccm}}, 
  \code{\link{dccm.xyz}}, \code{\link{dccm.nma}}, \code{\link{dccm.enma}}.
}
\examples{
\donttest{

##-- Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## Select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb, resno=c(24:27,85:90), elety='CA')

## lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

## Do PCA
pca <- pca.xyz(xyz)

## DCCM: only use first 10 PCs
cij <- dccm(pca, pc = c(1:10))

## Plot DCCM
plot(cij)

## DCCM: remove first 10 PCs
cij <- dccm(pca, pc = -c(1:10))

## Plot DCCM
plot(cij)

}
}
\keyword{ analysis }
