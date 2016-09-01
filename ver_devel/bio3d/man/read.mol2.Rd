\name{read.mol2}
\alias{read.mol2}
\alias{print.mol2}
\title{ Read MOL2 File }
\description{
  Read a Tripos MOL2 file
}
\usage{
read.mol2(file, maxlines = -1L)

\method{print}{mol2}(x, \dots)
}
\arguments{
  \item{file}{ a single element character vector containing the name of
    the MOL2 file to be read. }
  \item{maxlines}{ the maximum number of lines to read before giving up
    with large files. Default is all lines. }
  \item{x}{ an object as obtained from \code{read.mol2}. }
  \item{...}{ additional arguments to \sQuote{print}. }
}
\details{
  Basic functionality to parse a MOL2 file. The current version reads
  and stores \sQuote{@<TRIPOS>MOLECULE}, \sQuote{@<TRIPOS>ATOM},
  \sQuote{@<TRIPOS>BOND} and \sQuote{@<TRIPOS>SUBSTRUCTURE} records. 

  In the case of a multi-molecule MOL2 file, each molecule will be stored
  as an individual \sQuote{mol2} object in a list. Conversely, if the multi-molecule
  MOL2 file contains identical molecules in different conformations
  (typically from a docking run), then the output will be one object
  with an \code{atom} and \code{xyz} component (xyz in
  matrix representation; row-wise coordinates).

  See examples for further details.
}
\value{
  Returns a list of molecules containing the following components:
  \item{atom}{ a data frame containing all atomic coordinate ATOM data, 
    with a row per ATOM and a column per record type.  See below for details
    of the record type naming convention (useful for accessing
    columns). }
  \item{bond}{ a data frame containing all atomic bond information. }
  \item{substructure}{ a data frame containing all substructure information. }
  \item{xyz}{ a numeric matrix of ATOM coordinate data. }
  \item{info}{ a numeric vector of MOL2 info data. }
  \item{name}{ a single element character vector containing the molecule name. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
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
  Reisude number \dQuote{resno},
  Atom type \dQuote{elety},
  Residue name \dQuote{resid},
  Atom charge \dQuote{charge},
  Status bit \dQuote{statbit},

  For \code{bond} list components the column names are:
  Bond identifier \dQuote{id},
  number of the atom at one end of the bond\dQuote{origin},
  number of the atom at the other end of the bond \dQuote{target},
  the SYBYL bond type \dQuote{type}.

  For \code{substructure} list components the column names are:
  substructure identifier \dQuote{id},
  substructure name \dQuote{name},
  the ID number of the substructure's root atom \dQuote{root_atom},
  the substructure type \dQuote{subst_type},
  the type of dictionary associated with the substructure \dQuote{dict_type},
  the chain to which the substructre belongs \dQuote{chain},
  the subtype of the chain \dQuote{sub_type},
  the number of inter bonds \dQuote{inter_bonds},
  status bit \dQuote{status}. 
    
  See examples for further details.
}
\seealso{
  \code{\link{write.mol2}}, \code{\link{atom.select.mol2}}, 
  \code{\link{trim.mol2}}, \code{\link{as.pdb.mol2}}
  \code{\link{read.pdb}}
}
\examples{
cat("\n")
\dontrun{
## Read a single entry MOL2 file
## (returns a single object)
mol <- read.mol2( system.file("examples/aspirin.mol2", package="bio3d") )

## Short summary of the molecule
print(mol)

## ATOM records
mol$atom

## BOND records
mol$bond

## Print some coordinate data
head(mol$atom[, c("x","y","z")])

## Or coordinates as a numeric vector
#head(mol$xyz)

## Print atom charges
head(mol$atom[, "charge"])

## Convert to PDB
pdb <- as.pdb(mol)


## Read a multi-molecule MOL2 file
## (returns a list of objects)
#multi.mol <- read.mol2("zinc.mol2")

## Number of molecules described in file
#length(multi.mol)

## Access ATOM records for the first molecule
#multi.mol[[1]]$atom

## Or coordinates for the second molecule
#multi.mol[[2]]$xyz



## Process output from docking (e.g. DOCK)
## (typically one molecule with many conformations)
## (returns one object, but xyz in matrix format)
#dock <- read.mol2("dock.mol2")

## Reference PDB file (e.g. X-ray structure)
#pdb <- read.pdb("dock_ref.pdb")

## Calculate RMSD of docking modes
#sele <- atom.select(dock, "noh")
#rmsd(pdb$xyz, dock$xyz, b.inds=sele$xyz)
}
}
\keyword{ IO }
