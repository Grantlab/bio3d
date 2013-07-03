\name{torsion.pdb}
\alias{torsion.pdb}
\title{ Calculate Mainchain and Sidechain Torsion/Dihedral Angles }
\description{
  Calculate all torsion angles for a given protein PDB structure object.
}
\usage{
torsion.pdb(pdb)
}
\arguments{
  \item{pdb}{ a PDB structure object as obtained from
    function \code{read.pdb}. }
}
\details{
  The conformation of a polypeptide chain can be usefully described in
  terms of angles of internal rotation around its constituent bonds. See
  the related \code{torsion.xyz} function, which is called by this
  function, for details.
}
\value{
  Returns a list object with the following components:
  \item{phi}{ main chain torsion angle for atoms C,N,CA,C. }
  \item{psi}{ main chain torsion angle for atoms N,CA,C,N. }
  \item{omega}{ main chain torsion angle for atoms CA,C,N,CA. }
  \item{alpha}{ virtual torsion angle between consecutive C-alpha atoms. }
  \item{chi1}{ side chain torsion angle for atoms N,CA,CB,*G. }
  \item{chi2}{ side chain torsion angle for atoms CA,CB,*G,*D. }
  \item{chi3}{ side chain torsion angle for atoms CB,*G,*D,*E. }
  \item{chi4}{ side chain torsion angle for atoms *G,*D,*E,*Z. }
  \item{chi5}{ side chain torsion angle for atoms *D,*E,*Z, NH1. }
  \item{coords}{ numeric matrix of \sQuote{justified} coordinates. }
  \item{tbl}{ a numeric matrix of psi, phi and chi torsion angles. }
  
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ 
  For the protein backbone, or main-chain atoms, the partial double-bond
  character of the peptide bond between \sQuote{C=N} atoms severely restricts
  internal rotations.  In contrast, internal rotations around the single
  bonds between \sQuote{N-CA} and \sQuote{CA-C} are only restricted by
  potential steric collisions. Thus, to a good approximation, the
  backbone conformation of each residue in a given polypeptide chain can
  be characterised by the two angles phi and psi.

  Sidechain conformations can also be described by angles of internal
  rotation denoted chi1 up to chi5 moving out along the sidechain.
}
\seealso{ \code{\link{torsion.xyz}}, \code{\link{read.pdb}},
  \code{\link{dssp}}, \code{\link{stride}}. }
\examples{
##-- PDB torsion analysis
pdb <- read.pdb( "1bg2" )
tor <- torsion.pdb(pdb)
head(tor$tbl)

## basic Ramachandran plot
plot(tor$phi, tor$psi)

## torsion analysis of a single coordinate vector 
inds <- atom.select(pdb,"calpha")
tor.ca <- torsion.xyz(pdb$xyz[inds$xyz], atm.inc=1)

##-- Compare two PDBs to highlight interesting residues
aln <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))
m <- read.fasta.pdb(aln)
a <- torsion.xyz(m$xyz[1,],1)
b <- torsion.xyz(m$xyz[2,],1)
d <- wrap.tor(a-b)
plot(m$resno[1,],d, typ="h")
}
\keyword{ utilities }
