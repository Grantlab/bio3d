\name{torsion.xyz}
\alias{torsion.xyz}
\title{ Calculate Torsion/Dihedral Angles }
\description{
  Defined from the Cartesian coordinates of four successive atoms
  (A-B-C-D) the torsion or dihedral angle is calculated about an axis
  defined by the middle pair of atoms (B-C).
}
\usage{
torsion.xyz(xyz, atm.inc = 4)
}
\arguments{
  \item{xyz}{ a numeric vector of Cartisean coordinates. }
  \item{atm.inc}{ a numeric value indicating the number of atoms to
    increment by between successive torsion evaluations (see below). }
}
\details{
  The conformation of a polypeptide or nucleotide chain can be usefully
  described in terms of angles of internal rotation around its
  constituent bonds.

  If a system of four atoms A-B-C-D is projected onto a plane normal to
  bond B-C, the angle between the projection of A-B and the projection
  of C-D is described as the torsion angle of A and D about bond B-C.

  By convention angles are measured in the range -180 to +180, rather
  than from 0 to 360, with positive values defined to be in the
  clockwise direction.

  With \code{atm.inc=1}, torsion angles are calculated for each set of
  four successive atoms contained in \code{xyz} (i.e. moving along one
  atom, or three elements of \code{xyz}, between sucessive
  evaluations).  With \code{atm.inc=4}, torsion angles are calculated
  for each set of four successive non-overlapping atoms contained in
  \code{xyz} (i.e. moving along four atoms, or twelve elements of
  \code{xyz}, between sucessive evaluations).
  
}
\value{
  A numeric vector of torsion angles.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Karim ElSawy }
\note{ Contributions from Barry Grant. }
\seealso{ \code{\link{torsion.pdb}}, \code{\link{pca.tor}},
  \code{\link{wrap.tor}}, \code{\link{read.pdb}}, \code{\link{read.dcd}}.
}
\examples{
## Calculate torsions for cis & trans conformers
xyz <- rbind(c(0,-0.5,0,1,0,0,1,1,0,0,1.5,0),
             c(0,-0.5,0,1,0,0,1,1,0,2,1.5,0)-3)
cis.tor   <- torsion.xyz( xyz[1,] )
trans.tor <- torsion.xyz( xyz[2,] )
apply(xyz, 1, torsion.xyz)


plot(range(xyz), range(xyz), xlab="", ylab="", typ="n", axes=FALSE)
  apply(xyz, 1, function(x){
    lines(matrix(x, ncol=3, byrow=TRUE), lwd=4)
    points(matrix(x, ncol=3, byrow=TRUE), cex=2.5,
           bg="white", col="black", pch=21) } )
  
  text( t(apply(xyz, 1, function(x){
    apply(matrix(x, ncol=3, byrow=TRUE)[c(2,3),], 2, mean) })),
       labels=c(0,180), adj=-0.5, col="red")


\dontrun{
##-- PDB torsion analysis
pdb <- read.pdb("1bg2")
tor <- torsion.pdb(pdb)
## basic Ramachandran plot
plot(tor$phi, tor$psi)


## torsion analysis of a single coordinate vector 
inds <- atom.select(pdb,"calpha")
tor.ca <- torsion.xyz(pdb$xyz[inds$xyz], atm.inc=3)


##-- Compare two PDBs to highlight interesting residues
aln <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))
m <- read.fasta.pdb(aln)
a <- torsion.xyz(m$xyz[1,],1)
b <- torsion.xyz(m$xyz[2,],1)
## Note the periodicity of torsion angles
d <- wrap.tor(a-b)
plot(m$resno[1,],d, typ="h")
}

}
\keyword{ utilities }
