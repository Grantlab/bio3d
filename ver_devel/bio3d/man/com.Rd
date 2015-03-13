\name{com}
\alias{com}
\alias{com.pdb}
\alias{com.xyz}
\title{ Center of Mass }
\description{
  Calculate the center of mass of a PDB object.
}
\usage{
com(...)

\method{com}{pdb}(pdb, inds=NULL, use.mass=TRUE, ...)

\method{com}{xyz}(xyz, mass=NULL, ...)
}
\arguments{
  \item{pdb}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{pdb} upon which the calculation
    should be based.}
  \item{use.mass}{ logical, if TRUE the calculation will be mass
    weighted (center of mass). }
  \item{...}{ additional arguments to \code{atom2mass}. }
  \item{xyz}{ a numeric vector or matrix of Cartesian coordinates
    (e.g. an object of type \code{xyz}). }
  \item{mass}{ a numeric vector containing the masses of each atom in
    \code{xyz}. }
}
\details{
  This function calculates the center of mass of the provided PDB
  structure / Cartesian coordiantes. Atom names found in standard amino
  acids in the PDB are mapped to atom elements and their corresponding
  relative atomic masses.  

  In the case of an unknown atom name \code{elety.custom} and
  \code{mass.custom} can be used to map an atom to the correct
  atomic mass. See examples for more details.

  Alternatively, the atom name will be mapped automatically to the
  element corresponding to the first character of the atom name. Atom
  names starting with character \code{H} will be mapped to hydrogen
  atoms.
}
\value{
  Returns the Cartesian coordinates at the center of mass.
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{read.pdb}}, \code{\link{atom2mass}}
}
\examples{
\donttest{
# PDB server connection required - testing excluded

## Stucture of PKA:
pdb <- read.pdb("3dnd")

## Center of mass:
com(pdb)

## Center of mass of a selection
inds <- atom.select(pdb, chain="I")
com(pdb, inds)

## using XYZ Cartesian coordinates
xyz <- pdb$xyz[, inds$xyz]
com.xyz(xyz)

## with mass weighting
com.xyz(xyz, mass=atom2mass(pdb$atom[inds$atom, "elety"]) )

}

\dontrun{
## Unknown atom names
pdb <- read.pdb("3dnd")
inds <- atom.select(pdb, resid="LL2")
mycom <- com(pdb, inds, rescue=TRUE)
#warnings()


## Map atom names manually
pdb <- read.pdb("3RE0")
inds <- atom.select(pdb, resno=201)

myelety  <- data.frame(name = c("CL2","PT1","N1","N2"), symb = c("Cl","Pt","N","N"))
mymasses <- data.frame(symb = c("Cl","Pt"), mass = c(35.45, 195.08))
mycom    <- com(pdb, inds, elety.custom=myelety, mass.custom=mymasses)
}
}
\keyword{ utilities }
