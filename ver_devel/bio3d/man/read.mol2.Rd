\name{read.mol2}
\alias{read.mol2}
\title{ Read MOL2 File }
\description{
  Read a Sybyl MOL2 file
}
\usage{
read.mol2(file, maxlines = -1L)
}
\arguments{
  \item{file}{ a single element character vector containing the name of
    the MOL2 file to be read. }
  \item{maxlines}{ the maximum number of lines to read before giving up
    with large files. Default is all lines. }
}
\details{
  Basic functionality to parse a MOL2 file. The current version omits
  bond information, and only \code{@<TRIPOS>MOLECULE} and
  \code{@<TRIPOS>ATOM} records are stored.

  In the case of a multi-molecule MOL2 file, each molecule will be stored
  as an individual object in a list. Conversely, if the multi-molecule
  MOL2 file contains identical molecules in different conformations
  (typically a dockin run), then the output will be one object
  with an \code{atom} and \code{xyz} component (xyz in
  matrix representation; row-wise coordinates).

  See examples for further details.
}
\value{
  Returns a list of molecules containing the following components:
  \item{atom}{ a character matrix containing all atomic coordinate ATOM data, 
  with a row per ATOM and a column per record type.  See below for details
  of the record type naming convention (useful for accessing columns). }
  \item{xyz}{ a numeric vector or matrix of ATOM coordinate data. }
  \item{info}{ a numeric vector of MOL2 info data. }
  \item{name}{ a single element character vector containing the molecule name. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of the MOL2 format see:\cr
  \url{http://www.tripos.com/data/support/mol2.pdf}.
}
\author{ Lars Skjaerven }
\note{ 
  For \code{atom} list components the column names can be
  used as a convenient means of data access, namely:
  Atom serial number \dQuote{eleno},
  Atom name \dQuote{elena},
  Orthogonal coordinates \dQuote{x},
  Orthogonal coordinates \dQuote{y},
  Orthogonal coordinates \dQuote{z},
  Atom type \dQuote{elety},
  Residue name \dQuote{resid},
  Atom charge \dQuote{charge},
  Status bit \dQuote{statbit},
  See examples for further details.
}
\seealso{ \code{\link{atom.select}}, \code{\link{read.pdb}} }
\examples{
\dontrun{
## Read a single entry MOL2 file
## (returns a single object)
mol <- read.mol2("single.mol2")

## ATOM records
mol$atom

## Print some coordinate data
head(mol$atom[, c("x","y","z")])

## Or coordinates as a numeric vector
head(mol$xyz)

## Print atom charges
head(mol$atom[, "charge")])



## Read a multi-molecule MOL2 file
## (returns a list of objects)
multi.mol <- read.mol2("zinc.mol2")

## Number of molecules described in file
length(multi.mol)

## Access ATOM records for the first molecule
multi.mol[[1]]$atom

## Or coordinates for the second molecule
multi.mol[[2]]$xyz



## Process output from docking (e.g. DOCK)
## (typically one molecule with many conformations)
## (returns one object, but xyz in matrix format)
dock <- read.mol2("dock.mol2")

## Reference PDB file (e.g. X-ray structure)
pdb <- read.pdb("dock_ref.pdb", het2atom=T)

## Calculate RMSD of docking modes
sele <- atom.select(dock, "noh")
rmsd(pdb$xyz, dock$xyz, b.inds=sele$xyz)
}
}
\keyword{ IO }
