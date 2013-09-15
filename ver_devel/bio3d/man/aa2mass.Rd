\name{aa2mass}
\alias{aa2mass}
\title{ Amino Acid Residues to Mass Converter }
\description{
  Convert a sequence of amino acid residue names to mass.
}
\usage{
aa2mass(pdb, inds=NULL, mass.custom=NULL, addter=TRUE)
}
\arguments{
  \item{pdb}{ a character vector containing the atom names to convert
    to atomic masses. Alternatively, a object of type \code{pdb} can be
    provided. }
  \item{inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{pdb} upon which the calculation
    should be based. }
  \item{mass.custom}{ a list of amino acid residue names and their
    corresponding masses. }
  \item{addter}{ logical, if TRUE terminal atoms are added to final
    masses. }
 }
\details{
  This function converts amino acid residue names to their corresponding
  masses. In the case of a non-standard amino acid residue name
  \code{mass.custom} can be used to map the residue to the correct mass.

  See examples for more details.
}
\note{
  When object of type \code{pdb} is provided, non-calpha atom records
  are omitted from the selection.
}
\value{
  Returns a numeric vector of masses.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{atom.index}}, \code{\link{atom2mass}},
  \code{\link{aa.index}}}
\examples{
resi.names <- c("LYS", "ALA", "CYS", "HIS")
masses <- aa2mass(resi.names)

\dontrun{
## Fetch atomic masses in a PDB object
pdb <- read.pdb("3dnd")
masses <- aa2mass(pdb)

## or
masses <- atom2mass(pdb$atom[1:10,"resid"])


## Map residue names manually
pdb <- read.pdb("prot.amber.pdb")

mass.cust <- list("HIE"=137.141527, "HID"=137.141527)
masses <- aa2mass(pdb, mass.custom=mass.cust)
}
}
\keyword{ utilities }
