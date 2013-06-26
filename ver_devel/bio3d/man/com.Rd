\name{com}
\alias{com}
\alias{com.xyz}
\title{ Center of Mass }
\description{
  Calculate the center of mass of a PDB object.
}
\usage{
com(pdb, inds=NULL, use.mass=TRUE,
    mass.custom=NULL, elety.custom=NULL, rescue=FALSE)
com.xyz(xyz, mass=NULL)
}
\arguments{
  \item{pdb}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{pdb} upon which the calculation
    should be based.}
  \item{use.mass}{ logical, if TRUE the calculation will be mass
    weighted (center of mass). }
  \item{mass.custom}{ a list of atom elements and their corresponding masses. }
  \item{elety.custom}{ a list of atom names to element symbol mapping. }
  \item{rescue}{ logical, if TRUE the atom element will be mapped to the
    first character of the atom name. }
  \item{xyz}{ a numeric vector of Cartesian coordinates. }
  \item{mass}{ a numeric vector containing the masses of each atom in
    \code{xyz}. }
}
\details{
  This function calculates the center of mass of the provided PDB
  structure. Atom names found in standard amino acids in the PDB are
  mapped to atom elements and their corresponding relative atomic masses. 

  In the case of an unknown atom name \code{elety.custom} and
  \code{mass.custom} can be used to map an atom to the correct
  atomic mass. See examples for more details.

  Alternatively, the atom name will be mapped automatically to the
  element corresponding to the first character of the atom name. Atom
  names starting with character \code{H} will be mapped to hydrogen
  atoms. 
 
  The function reports the atomic masses in which the calculation is
  based and the corresponding center of mass. 
}
\value{
  Returns a list with the following components:
  \item{mass}{ a numeric vector of atomic masses. }
  \item{xyz}{ Cartesian coordinates of the center of mass. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjærven }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}} }
\examples{

## Stucture of PKA:
pdb <- read.pdb("3dnd")

## Center of mass:
mycom <- com(pdb)
write.pdb(pdb = NULL, xyz=mycom$xyz, file="dummy-1.pdb", append=F)

## Center of mass of a selection
inds <- atom.select(pdb, chain="I")
mycom <- com(pdb, inds)

## Unknown atom names
pdb <- read.pdb("3dnd", het2atom=TRUE)
inds <- atom.select(pdb, resid="LL2")
mycom <- com(pdb, inds, rescue=TRUE)
warnings()


## Map atom names manually
pdb <- read.pdb("3RE0", het2atom=T)
inds <- atom.select(pdb, resno=201)

elety.cust <- list("CL2"="Cl", "PT1"="Pt")
mass.cust <- list("Cl"=35.45, "Pt"=195.08)

mycom <- com(pdb, inds, elety.custom=elety.cust, mass.custom=mass.cust,
             rescue=TRUE)

## Load data for atomic mass and atom names to element mapping
data(atom.index)

}
\keyword{ utilities }