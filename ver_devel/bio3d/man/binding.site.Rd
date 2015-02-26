\name{binding.site}
\alias{binding.site}
\title{ Binding Site Residues }
\description{
  Determines the interacting residues between two PDB entities.
}
\usage{
binding.site(a, b=NULL, a.inds=NULL, b.inds=NULL, cutoff=5,
            hydrogens=TRUE, byres=TRUE, verbose=FALSE)
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
  \item{cutoff}{ distance cutoff }
  \item{hydrogens}{ logical, if FALSE hydrogen atoms are omitted from
    the calculation. }
  \item{byres}{ logical, if TRUE all atoms in a contacting residue is
    returned. }
  \item{verbose}{ logical, if TRUE details of the selection are printed. }
}
\details{
  This function reports the residues of \code{a} closer than a cutoff to
  \code{b}. This is a wrapper function calling the underlying function
  \code{dist.xyz}.

  If \code{b=NULL} then \code{b.inds} should be elements of \code{a}
  upon which the calculation is based (typically chain A and B of the
  same PDB file).

  If \code{b=a.inds=b.inds=NULL} the function will use
  \code{\link{atom.select}} with arguments \code{"protein"} and
  \code{"ligand"} to determine receptor and ligand, respectively.
}
\value{
  Returns a list with the following components:
  \item{inds}{ object of class \code{select} with \code{atom} and
  \code{xyz} components. }
  \item{inds$atom}{ atom indices of \code{a}. }
  \item{inds$xyz}{ xyz indices of \code{a}. }
  \item{resnames}{ a character vector of interacting residues. }
  \item{resno}{ a numeric vector of interacting residues numbers. }
  \item{chain}{ a character vector of the associated chain identifiers
    of \code{"resno"}. }
  \item{call}{ the matched call. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}}, \code{\link{dm}} }
\examples{
   \donttest{
    # PDB server connection required - testing excluded

     pdb <- read.pdb('3dnd')

     ## automatically identify 'protein' and 'ligand'
     bs <- binding.site(pdb)

     bs$resnames
     #pdb$atom[bs$inds$atom, ]

     # provide indices
     rec.inds <- atom.select(pdb, chain="A", resno=1:350)
     lig.inds <- atom.select(pdb, chain="A", resno=351)
     bs <- binding.site(pdb, a.inds=rec.inds, b.inds=lig.inds)
   }


   \dontrun{  
     # Interaction between peptide and protein
     rec.inds <- atom.select(pdb, chain='A', resno=c(1:350))
     lig.inds <- atom.select(pdb, chain='I', resno=c(5:24))
     bs <- binding.site(pdb, a.inds=rec.inds, b.inds=lig.inds)
   }

   \donttest{
    # Redundant testing excluded

     # Interaction between two PDB entities   
     #rec <- read.pdb("receptor.pdb")
     #lig <- read.pdb("ligand.pdb")
     rec <- trim.pdb(pdb, inds=rec.inds)
     lig <- trim.pdb(pdb, inds=lig.inds)
     bs <- binding.site(rec, lig, hydrogens=FALSE)
   }
}
\keyword{ utilities }
