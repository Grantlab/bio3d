\name{dccm}
\alias{dccm}
\title{ DCCM: Dynamical Cross-Correlation Matrix }
\description{
  Determine the cross-correlations of atomic displacements.
}
\usage{
dccm(xyz, reference = apply(xyz, 2, mean), ncore=1, nseg.scale=1)
}
\arguments{
  \item{xyz}{ a numeric matrix of Cartesian coordinates with a row per
    structure/frame. }
  \item{reference}{ The reference structure about which displacements are
    analysed. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. } 
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
  The extent to which the atomic fluctuations/displacements of a system are
  correlated with one another can be assessed by examining the magnitude
  of all pairwise cross-correlation coefficients (see McCammon and Harvey,
  1986).

  This function returns a matrix of all atom-wise cross-correlations
  whose elements, Cij, may be displayed in a graphical representation
  frequently termed a dynamical cross-correlation map, or DCCM.

  If Cij = 1 the fluctuations of atoms i and j are completely correlated
  (same period and same phase), if Cij = -1 the fluctuations of atoms i
  and j are completely anticorrelated (same period and opposite phase),
  and if Cij = 0 the fluctuations of i and j are not correlated.

  Typical characteristics of DCCMs include a line of strong
  cross-correlation along the diagonal, cross-correlations emanating
  from the diagonal, and off-diagonal cross-correlations. The high
  diagonal values occur where i = j, where Cij is always equal to
  1.00. Positive correlations emanating from the diagonal indicate
  correlations between contiguous residues, typically within a secondary
  structure element or other tightly packed unit of structure.
  Typical secondary structure patterns include a triangular pattern for
  helices and a plume for strands.  Off-diagonal positive and negative
  correlations may indicate potentially interesting correlations between
  domains of non-contiguous residues.
}
\value{
  Returns a cross-correlation matrix.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  McCammon, A. J. and Harvey, S. C. (1986) \emph{Dynamics of
  Proteins and Nucleic Acids}, Cambridge University Press, Cambridge.
}
\author{ Gisle Saelensminde }
\note{
  This function is currently very basic i.e. inefficient and \bold{SLOW}. 
}
\seealso{ \code{\link{cor}} for examining xyz cross-correlations,
  \code{\link{pca.xyz}}. }
\examples{

\dontrun{
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

## DCCM (slow to run so restrict to Calpha)
cij <- dccm(xyz)

## Plot DCCM
plot(cij)

## Or
library(lattice)
contourplot(cij, region = TRUE, labels=F, col="gray40",
            at=c(-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1),
            xlab="Residue No.", ylab="Residue No.",
            main="DCCM: dynamic cross-correlation map")
}

}
\keyword{ utilities }
