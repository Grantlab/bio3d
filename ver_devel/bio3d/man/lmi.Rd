\name{lmi}
\alias{lmi}
\title{ LMI: Linear Mutual Information Matrix }
\description{
  Calculate the linear mutural information correlations of atomic displacements.
}
\usage{
lmi(trj, grpby=NULL, ncore=1)
}
\arguments{
  \item{trj}{ a numeric matrix of Cartesian coordinates with a row per structure/frame. }
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of \code{trj} that should be considered as a group
    (e.g. atoms from a particular residue). }
  \item{ncore}{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
}
\details{
  The correlation of the atomic fluctuations of a system can be assessed by the Linear 
  Mutual Information (LMI) and the LMI has no unwanted dependency on the relative orientation 
  of the fluctuations which the Pearson coefficient suffers from.

  This function returns a matrix of all atom-wise linear mutual information whose elements are
  denoted as Cij. If Cij = 1, the fluctuations of atoms i and j are completely correlated and if
  Cij = 0, the fluctuations of atoms i and j are not correlated.
}
\value{
  Returns a linear mutual information matrix.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  Lange, O.F. and Grubmuller, H. (2006) \emph{PROTEINS: Structure, Function, and Bioinformatics} \bold{62}:1053--1061.
}
\author{ Hongyang Li & Barry Grant}
\examples{
##-- Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb,"///24:27,85:90///CA/")

## lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

## LMI matrix (slow to run so restrict to Calpha)
cij <- lmi(xyz)

## Plot LMI matrix
#plot(cij)
col.scale <- colorRampPalette(c("gray95", "cyan"))(5)
plot(cij, at=seq(0.4,1, length=5), col.regions=col.scale)

}


