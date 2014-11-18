\name{atom.select}
\alias{atom.select}
\alias{atom.select.pdb}
\alias{atom.select.prmtop}
\alias{print.select}
\title{ Atom Selection from PDB and PRMTOP Structure Objects }
\description{
  Return the \sQuote{atom} and \sQuote{xyz} coordinate indices of
  \sQuote{pdb} or \sQuote{prmtop} structure objects corresponding to the
  intersection of a hierarchical selection. 
}
\usage{
atom.select(\dots)

\method{atom.select}{pdb}(pdb, string=NULL, chain=NULL, resno=NULL, resid=NULL,
            eleno=NULL, elety=NULL, type=NULL, verbose=TRUE, rm.insert=FALSE, \dots)

\method{atom.select}{prmtop}(prmtop, ...)

\method{print}{select}(x, \dots)
}
\arguments{
  \item{\dots}{ arguments passed to \code{atom.select.pdb},
    \code{atom.select.prmtop}, or \code{print}. }
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{string}{ a character selection string with the following syntax:\cr
    \code{/segid/chain/resno/resid/eleno/elety/}. Or a single selection keyword from 
    \code{calpha} \code{cbeta} \code{backbone} \code{protein} \code{notprotein} \code{ligand} \code{water} \code{notwater} \code{h} \code{noh} }
  \item{chain}{ a character vector of chain identifiers. }
  \item{resno}{ a numeric or character vector of residue numbers. }
  \item{resid}{ a character vector of residue name identifiers. }
  \item{eleno}{ a numeric or character vector of element numbers. }
  \item{elety}{ a character vector of atom names. }
  \item{type}{ a single element character vector for selecting \sQuote{ATOM}
     or \sQuote{HETATM} record types. }
  \item{verbose}{ logical, if TRUE details of the selection are printed. }
  \item{rm.insert}{ logical, if TRUE insert ATOM records from the
    \code{pdb} object are ignored. }
  \item{prmtop}{ a structure object of class \code{"prmtop"}, obtained from
    \code{\link{read.prmtop}}. }
  \item{x}{ a atom.select object as obtained from
    \code{\link{atom.select}}.  }
}
\details{
  This function allows for the selection of atom and coordinate data
  corresponding to the intersection of various input criteria. 

  Input selection criteria include selection \code{string} keywords (such as 
  \code{"calpha"}, \code{"backbone"}, \code{"protein"}, \code{"ligand"}, etc.)
  and individual named selection components (including \sQuote{chain}, 
  \sQuote{resno}, \sQuote{resid}, \sQuote{elety} etc.). For example,
  \code{atom.select(pdb, "calpha")} or \code{atom.select(pdb,
  elety="CA")}, which in this case are equivalent.

  Alternatively, a single element character vector containing a hierarchical
  \sQuote{selection string}, \code{string}, with a strict format composed of 
  six sections separated by a \sQuote{/} character can be used.
  
  Each section of this \sQuote{selection string} corresponds to a
  different level in the PDB structure hierarchy, namely:
  (1) segment identifier,
  (2) chain identifier,
  (3) residue number,
  (4) residue name,
  (5) element number, and
  (6) element name.

  For example, the string \code{//A/65:143///CA/} selects all C-alpha
  atoms from residue numbers 65 to 143, of chain A.

  A simpler alternative would be \code{atom.select(pdb, chain="A", resno=65:143,
    elety="CA")}. 

  More typical use is through keyword \code{string} shortcuts
  \code{"calpha"}, \code{"back"}, \code{"backbone"}, \code{"cbeta"},
  \code{"protein"}, \code{"notprotein"}, \code{"ligand"}, \code{"water"}, 
  \code{"notwater"}, \code{"h"} and \code{"noh"}. 

  Note that keyword \code{string} shortcuts can be combined with individual
  selection components, e.g. \code{atom.select(pdb, "protein",
  chain="A")}. See below for further examples.  

  In addition, the \code{combine.sel} function can further combine atom
  selections using \sQuote{and}, \sQuote{or}, or \sQuote{not} logical
  operations.   

  When called without selection criteria, \code{atom.select} will print a
  summary of \code{pdb} makeup.

  See also \code{atom.select.pdb} or \code{atom.select.prmtop} for more
  information. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{
  \code{\link{combine.sel}}, \code{\link{read.pdb}}, 
  \code{\link{write.pdb}}, \code{\link{read.prmtop}},
  \code{\link{read.crd}}, \code{\link{read.dcd}},
  \code{\link{read.ncdf}}.
}
\examples{

##- PDB example
# Read a PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

# Print structure summary
atom.select(pdb)

# Select all C-alpha atoms with residues numbers between 43 and 54
ca.inds <- atom.select(pdb, "calpha", resno=43:54)
print( pdb$atom[ ca.inds$atom, "resid" ] )
print( pdb$xyz[ ca.inds$xyz ] )

\dontrun{
# Select all C-alphas
ca.inds <- atom.select(pdb, "calpha")

# String examples (see above).
ca.inds <- atom.select(pdb, "///43:54///CA/")
inds <- atom.select(pdb, "///107:118///N,CA,C,O/")


##- PRMTOP example
prmtop <- read.prmtop("prot_solvated.prmtop")

## Atom selection 
ca.inds <- atom.select(prmtop, "calpha")

}
}
\keyword{utilities}

