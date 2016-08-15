\name{dccm.local.xyz}
\alias{dccm.local.xyz}
\title{
  Dynamical Cross-Correlation Matrix on local protein segments.
}
\description{
 dccm.local.xyz calculates the residue pairwise correlations using the
 RMSD values from a local fitting operation.
}
\usage{
dccm.local.xyz(trj, pdb, radius=9, basin.grp=NULL, write.average=FALSE, ncore=1)
}
\arguments{
  \item{trj}{ A numeric matrix of Cartesian coordinates with a row per
    structure/frame. }
  \item{pdb}{ A pdb object such as obtained by read.pdb(). }
  \item{radius}{ Numeric. Corresponds to the radius value (in A) for the
    sphere used in the neighboring residues determination. }
  \item{basin.grp}{ Numeric. A vector matching the number of frames in
    the trj object assigning the membership of each frame for the
    average structure calculation. }
  \item{write.average}{ Logic. If TRUE, average structures are written
    in pdb format. }
  \item{ncore}{ Numeric. Indicates the number of core to use. }
}
\value{
  A list containing two matrices, $rmsd.matrix corresponds to the
  residue RMSD matrix, $correlation.matrix to the pairwise residue correlations matrix.
}
\details{
  First, for each frame a reference structure is defined using either the
  previous trajectory frame (default) or an average structure from trajectory segments (through
 the frame membership vector frame.grp). Then, for each residue in each
 frame, the reference structure is used for the determination of the
 neighboring residues included in a sphere defined by the input
 radius. The neighbors of each residue are then used to perform the local fit
 and to calculate the local RMSD. Finally, on the colums of the obtained RMSD
 matrix (whose dimensions are (M-1)xN by default or MxN if the basin.grp
 vector is specified, where M is the number of frames
 and N the number of residues) the pearson correlations are calculated,
 resulting in a NxN correlation matrix. For more details, see the
 reference Balabin, I.A. et al. (2009).
  
}  

\references{
  Balabin, I.A. et al. (2009) \emph{PNAS} \bold{106(34)}, 14253--14258.  

  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}

\author{ Guido Scarabelli and Xin-Qiu Yao }

\seealso{
 \code{\link{cmap}}, \code{\link{fit.xyz}}, \code{\link{rmsd}}, \code{\link{cor}}, \code{\link{dccm.xyz}}.
}


\examples{
##-- Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

local.correlations <- dccm.local.xyz(trj, pdb, ncore=1)

plot.dccm(local.correlations$correlation.matrix)

}
