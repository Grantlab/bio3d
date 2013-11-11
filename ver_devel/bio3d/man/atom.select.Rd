\name{atom.select}
\alias{atom.select}
\title{ Atom Selection From PDB Structure }
\description{
  Return the atom and xyz coordinates indices of a \sQuote{pdb} structure
  object corresponding to the intersection of a hierarchical selection.
}
\usage{
atom.select(pdb, string=NULL, chain=NULL, resno=NULL, resid=NULL,
eleno=NULL, elety=NULL, verbose=TRUE, rm.insert=FALSE)
}
\arguments{
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
  \item{verbose}{ logical, if TRUE details of the selection are printed. }
  \item{rm.insert}{ logical, if TRUE insert ATOM records from the
    \code{pdb} object are ignored. }
}
\details{
  This function allows for the selection of atom and coordinate data
  corresponding to the intersection of input criteria (such as a
  \sQuote{chain}, \sQuote{resno}, \sQuote{elety} etc.  Or a hierarchical
  \sQuote{selection string} \code{string} with a strict format.
  
  The selection \code{string} should be a single element character
  vector containing a string composed of six sections separated by a
  \sQuote{/} character.
  
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

  A simpler alternative would be \code{chain="A", resno=65:143,
    elety="CA"}. In addition, the character \code{string} shortcuts
  \code{"calpha"}, \code{"back"}, \code{"backbone"}, \code{"cbeta"},
  \code{"h"} and \code{"noh"} may also be used. See below for examples.

  When called without a selection string, \code{atom.select} will print a
  summary of \code{pdb} makeup.
}
\value{
  Returns a list of class \code{"select"} with components:
  \item{atom }{atom indices of selected atoms.}
  \item{xyz }{xyz indices of selected atoms.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.pdb}}, \code{\link{write.pdb}}, \code{\link{read.dcd}} }
\examples{

# Read a PDB file
pdb <- read.pdb( "http://www.rcsb.org/pdb/files/1BG2.pdb" )

# Print structure summary
atom.select(pdb)

# Select all C-alpha atoms with residues numbers between 65 and 143
ca.inds <- atom.select(pdb, resno=65:143, elety="CA")
print( pdb$atom[ ca.inds$atom, "resid" ] )
print( pdb$xyz[ ca.inds$xyz ] )

\dontrun{
# Select all C-alphas
ca.inds <- atom.select(pdb, "calpha")

# String examples (see above).
ca.inds <- atom.select(pdb, "///65:143///CA/")
inds <- atom.select(pdb, "///130:142///N,CA,C,O/")
}
}
\keyword{utilities}

