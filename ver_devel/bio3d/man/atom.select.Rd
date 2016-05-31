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

\method{atom.select}{pdb}(pdb, string = NULL,
                          type  = NULL, eleno = NULL, elety = NULL,
                          resid = NULL, chain = NULL, resno = NULL,
                          segid = NULL, operator = "AND",
                          inverse = FALSE, value = FALSE, verbose=FALSE, \dots)

\method{atom.select}{prmtop}(prmtop, ...)

\method{print}{select}(x, \dots)
}
\arguments{
  \item{\dots}{ arguments passed to \code{atom.select.pdb},
    \code{atom.select.prmtop}, or \code{print}. }
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{string}{ a single selection keyword from \code{calpha}
    \code{cbeta} \code{backbone} \code{sidechain} \code{protein}
    \code{nucleic} \code{ligand} \code{water} \code{h} or \code{noh}. }
  \item{type}{ a single element character vector for selecting \sQuote{ATOM}
    or \sQuote{HETATM} record types. }
  \item{eleno}{ a numeric vector of element numbers. }
  \item{elety}{ a character vector of atom names. }
  \item{resid}{ a character vector of residue name identifiers. }
  \item{chain}{ a character vector of chain identifiers. }
  \item{resno}{ a numeric vector of residue numbers. }
  \item{segid}{ a character vector of segment identifiers. }
  \item{operator}{ a single element character specifying either the AND
    or OR operator by which individual selection components should be
    combined. Allowed values are \sQuote{"AND"} and \sQuote{"OR"}. }
  \item{verbose}{ logical, if TRUE details of the selection are printed. }
  \item{inverse}{ logical, if TRUE the inversed selection is retured
    (i.e. all atoms NOT in the selection). }
  \item{value}{ logical, if FALSE, vectors containing the (integer) indices of
    the matches determined by \code{atom.select} are returned, and if TRUE,
    a \code{pdb} object containing the matching atoms themselves is
    returned. }
  \item{prmtop}{ a structure object of class \code{"prmtop"}, obtained from
    \code{\link{read.prmtop}}. }
  \item{x}{ a atom.select object as obtained from
    \code{\link{atom.select}}.  }
}
\details{
  This function allows for the selection of atom and coordinate data
  corresponding to the intersection of various input criteria.

  Input selection criteria include selection \code{string} keywords (such as
  \code{"calpha"}, \code{"backbone"}, \code{"sidechain"},  \code{"protein"},
  \code{"nucleic"}, \code{"ligand"}, etc.) and individual named
  selection components (including \sQuote{chain},  \sQuote{resno},
  \sQuote{resid}, \sQuote{elety} etc.).

  For example, \code{atom.select(pdb, "calpha")} will return indices for
  all C-alpha (CA) atoms found in protein residues in the \code{pdb}
  object, \code{atom.select(pdb, "backbone")} will return indices
  for all protein N,CA,C,O atoms, and \code{atom.select(pdb, "cbeta")}
  for all protein N,CA,C,O,CB atoms. 

  Note that keyword \code{string} shortcuts can be combined with individual
  selection components, e.g. \code{atom.select(pdb, "protein", chain="A")}
  will select all protein atoms found in chain A.

  Selection criteria are combined according to the provided
  \code{operator} argument. The default operator \code{AND} (or \code{&}) will
  combine by intersection while \code{OR} (or \code{|}) will take the union.

  For example, \code{atom.select(pdb, "protein", elety=c("N", "CA", "C"),
  resno=65:103)} will select the N, CA, C atoms in the protein residues
  65 through 103, while \code{atom.select(pdb, "protein", resid="ATP",
  operator="OR")} will select all protein atoms as well as any ATP
  residue(s).

  Other \code{string} shortcuts include:
  \code{"calpha"}, \code{"back"}, \code{"backbone"}, \code{"cbeta"},
  \code{"protein"}, \code{"notprotein"}, \code{"ligand"},
  \code{"water"}, \code{"notwater"}, \code{"h"}, \code{"noh"},
  \code{"nucleic"}, and \code{"notnucleic"}.

  In addition, the \code{\link{combine.select}} function can further combine atom
  selections using \sQuote{AND}, \sQuote{OR}, or \sQuote{NOT} logical
  operations.

}
\note{
  Protein atoms are defined as any atom in a residue matching the
  residue name in the attached \code{aa.table} data frame. See
  \code{aa.table$aa3} for a complete list of residue names.

  Nucleic atoms are defined as all atoms found in residues with names A,
  U, G, C, T, I, DA, DU, DG, DC, DT, or DI.

  Water atoms/residues are defined as those with residue names H2O, OH2,
  HOH, HHO, OHH, SOL, WAT, TIP, TIP, TIP3, or TIP4.
}
\value{
  Returns a list of class \code{"select"} with the following components:
  \item{atom}{ a numeric matrix of atomic indices. }
  \item{xyz }{ a numeric matrix of xyz indices. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{
  \code{\link{read.pdb}}, \code{\link{as.select}},
  \code{\link{combine.select}}, \code{\link{trim.pdb}},
  \code{\link{write.pdb}}, \code{\link{read.prmtop}},
  \code{\link{read.crd}}, \code{\link{read.dcd}},
  \code{\link{read.ncdf}}.
}
\examples{

##- PDB example
# Read a PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

# Select protein atoms of chain A
atom.select(pdb, "protein", chain="A")

# Select all atoms except from the protein
atom.select(pdb, "protein", inverse=TRUE, verbose=TRUE)

# Select all C-alpha atoms with residues numbers between 43 and 54
sele <- atom.select(pdb, "calpha", resno=43:54, verbose=TRUE)

# Access the PDB data with the selection indices
print( pdb$atom[ sele$atom, "resid" ] )
print( pdb$xyz[ sele$xyz ] )

# Trim PDB to selection
ca.pdb <- trim.pdb(pdb, sele)

\dontrun{

##- PRMTOP example
prmtop <- read.prmtop("prot_solvated.prmtop")

## Atom selection
ca.inds <- atom.select(prmtop, "calpha")

}
}
\keyword{utilities}

