\name{aa2mass}
\alias{aa2mass}
\title{ Amino Acid Residues to Mass Converter }
\description{
  Convert a sequence of amino acid residue names to mass.
}
\usage{
aa2mass(pdb, inds=NULL, mass.custom=NULL, addter=TRUE, mmtk=FALSE)
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
  \item{mmtk}{ logical, if TRUE use the exact aminoacid residue masses as
    provided with the MMTK database (for testing purposes). }
}
\details{
  This function converts amino acid residue names to their corresponding
  masses. In the case of a non-standard amino acid residue name
  \code{mass.custom} can be used to map the residue to the correct
  mass. User-defined amino acid masses (with argument \code{mass.custom})
  will override mass entries obtained from the database. 

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
masses <- aa2mass(resi.names, addter=FALSE)

\dontrun{
## Fetch atomic masses in a PDB object
pdb <- read.pdb("3dnd")
masses <- aa2mass(pdb)

## or
masses <- atom2mass(pdb$atom[1:10,"resid"])

## Dealing with unconventional residues
pdb <- read.pdb("1xj0", het2atom=TRUE)

mass.cust <- list("CSX"=122.166)
masses <- aa2mass(pdb, mass.custom=mass.cust)
}
}
\keyword{ utilities }
