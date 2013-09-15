\name{atom2mass}
\alias{atom2mass}
\alias{atom2ele}
\title{ Atom Names to Mass Converter }
\description{
  Convert atom names to atomic mass.
}
\usage{
atom2mass(pdb, inds=NULL, mass.custom=NULL, elety.custom=NULL,
          grpby=NULL, rescue=TRUE)
atom2ele(pdb, inds=NULL, elety.custom=NULL, rescue=TRUE)
}
\arguments{
  \item{pdb}{ a character vector containing the atom names to convert
    to atomic masses. Alternatively, an object of type \code{pdb} can be
    provided. }
  \item{inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{pdb} upon which the calculation
    should be based (in effect only when a \code{pdb} object is provided). }
  \item{mass.custom}{ a list of atom elements and their corresponding masses. }
  \item{elety.custom}{ a list of atom names to element symbol mapping. }
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of \code{atoms} that should be considered as a group
    (e.g. atoms from a particular residue). }
  \item{rescue}{ logical, if TRUE the atom element will be mapped to the
    first character of the atom name. }
}
\details{
  This function converts atom names to their corresponding relative
  atomic masses. Atom names found in standard amino acids in the PDB are
  mapped to atom elements (with \code{atom2ele}) which finally are
  mapped to mass.

  In the case of an unknown atom name \code{elety.custom} and
  \code{mass.custom} can be used to map an atom to the correct
  atomic mass.

  Alternatively, the atom name will be mapped automatically to the
  element corresponding to the first character of the atom name. Atom
  names starting with character \code{H} will be mapped to hydrogen
  atoms.

  See examples for more details.
}
\value{
  Returns a numeric vector of atomic masses.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{atom.index}}, \code{\link{com}} }
\examples{
atom.names <- c("CA", "O", "N", "OXT")
masses <- atom2mass(atom.names)

\dontrun{
## Fetch atomic masses in a PDB object
pdb <- read.pdb("1hel")
masses <- atom2mass(pdb)

## or
masses <- atom2mass(pdb$atom[1:10,"elety"])

## Group and sum per residue
masses <- atom2mass(pdb, grpby=pdb$atom[,"resno"])


## Map atom names manually
pdb <- read.pdb("3RE0", het2atom=TRUE)
inds <- atom.select(pdb, resno=201)

elety.cust <- list("CL2"="Cl", "PT1"="Pt")
mass.cust <- list("Cl"=35.45, "Pt"=195.08)

masses <- atom2mass(pdb, inds, mass.custom=mass.cust, elety.custom=elety.cust)
}
}
\keyword{ utilities }
