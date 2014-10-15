\name{binding.site}
\alias{binding.site}
\title{ Binding Site Residues }
\description{
  Determines the interacting residues between two PDB entities.
}
\usage{
binding.site(a, b=NULL, a.inds=NULL, b.inds=NULL, cut=5, hydrogens=TRUE)
}
\arguments{
  \item{a}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{b}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. }
  \item{a.inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{a} upon which the calculation
    should be based.}
  \item{b.inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{b} upon which the calculation
    should be based.}
  \item{cut}{ distance cutoff }
  \item{hydrogens}{ logical, if FALSE hydrogen atoms are omitted from
    the calculation. }
}
\details{
  This function reports the residues of \code{a} closer than a cutoff to
  \code{b}. This is a wrapper function calling the underlying function
  \code{dm.xyz}.

  If \code{b=NULL} then \code{b.inds} should be elements of \code{a}
  upon which the calculation is based (typically chain A and B of the
  same PDB file).
}
\value{
  Returns a list with the following components:
  \item{atom.inds}{ atom indices of \code{a}. }
  \item{xyz.inds}{ xyz indices of \code{a}. }
  \item{resnames}{ a character vector of interacting residues. }
  \item{resno}{ a numeric vector of interacting residues numbers. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}}, \code{\link{dm}} }
\examples{
     \donttest{
     pdb <- read.pdb('3dnd')

     # Binding site residues
     rec.inds <- atom.select(pdb, string='//A/1:350////')
     lig.inds <- atom.select(pdb, string='//A/351////')
     bs <- binding.site(pdb, a.inds=rec.inds, b.inds=lig.inds)

     # Interaction between peptide and protein
     rec.inds <- atom.select(pdb, string='//A/1:350////')
     lig.inds <- atom.select(pdb, string='//I/5:24////')
     bs <- binding.site(pdb, a.inds=rec.inds, b.inds=lig.inds)
     }

     \dontrun{  
     # Interaction between two PDB entities
#     rec <- read.pdb("receptor.pdb")
#     lig <- read.pdb("ligand.pdb")
     rec <- trim.pdb(pdb, inds=rec.inds)
     lig <- trim.pdb(pdb, inds=lig.inds)
     bs <- binding.site(rec, lig, hydrogens=FALSE)
     }
}
\keyword{ utilities }
