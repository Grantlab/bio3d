\name{read.pqr}
\alias{read.pqr}
\title{ Read PQR File }
\description{
  Read a PQR coordinate file.
}
\usage{
read.pqr(file, maxlines = -1, multi = FALSE, rm.insert = FALSE,
         rm.alt = TRUE, verbose = TRUE)
}
\arguments{
  \item{file}{ the name of the PQR file to be read. }
  \item{maxlines}{ the maximum number of lines to read before giving up
    with large files. By default if will read up to the end of input on
    the connection. }
  \item{multi}{ logical, if TRUE multiple ATOM records are read for all
    models in multi-model files. }
  \item{rm.insert}{ logical, if TRUE PDB insert records are ignored. }
  \item{rm.alt}{ logical, if TRUE PDB alternate records are ignored. }
  \item{verbose}{ print details of the reading process. }
}
\details{
  PQR file format is basically the same as PDB format except for the fields of
  \code{o} and \code{b}. In PDB, these two fields are filled with \sQuote{Occupancy}
  and \sQuote{B-factor} values, respectively, with each field 6-column long.
  In PQR, they are atomic \sQuote{partial charge} and \sQuote{radii}
  values, respectively, with each field 8-column long.

  \code{maxlines} may require increasing for some large multi-model files.
  The preferred means of reading such data is via binary DCD format 
  trajectory files (see the \code{\link{read.dcd}} function).
}
\value{
  Returns a list of class \code{"pdb"} with the following components:
  \item{atom}{ a data.frame containing all atomic coordinate ATOM and
    HETATM data, with a row per ATOM/HETATM and a column per record type.
    See below for details of the record type naming convention (useful for
    accessing columns). }
  \item{helix }{ \sQuote{start}, \sQuote{end} and \sQuote{length} of H
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{sheet }{ \sQuote{start}, \sQuote{end} and \sQuote{length} of E
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{seqres }{ sequence from SEQRES field. }
  \item{xyz }{ a numeric matrix of class \code{"xyz"} containing the
    ATOM and HETATM coordinate data. }
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of PDB format (version3.3) see:\cr
  \url{http://www.wwpdb.org/documentation/format33/v3.3.html}.
}
\author{ Barry Grant }
\note{ 
  For both \code{atom} and \code{het} list components the column names can be
  used as a convenient means of data access, namely:
  Atom serial number \dQuote{eleno} ,
  Atom type \dQuote{elety},
  Alternate location indicator \dQuote{alt},
  Residue name \dQuote{resid},
  Chain identifier \dQuote{chain},
  Residue sequence number \dQuote{resno},
  Code for insertion of residues \dQuote{insert},
  Orthogonal coordinates \dQuote{x},
  Orthogonal coordinates \dQuote{y},
  Orthogonal coordinates \dQuote{z},
  Occupancy \dQuote{o}, and
  Temperature factor \dQuote{b}.
  See examples for further details.
}
\seealso{ \code{\link{atom.select}}, \code{\link{write.pqr}}, 
  \code{\link{read.pdb}}, \code{\link{write.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.fasta.pdb}},
  \code{\link{read.fasta}} }
\examples{
\donttest{
# PDB server connection required - testing excluded

# Read a PDB file and write it as a PQR file
pdb <- read.pdb( "4q21" )
outfile = file.path(tempdir(), "eg.pqr")
write.pqr(pdb=pdb, file = outfile)

# Read the PQR file
pqr <- read.pqr(outfile)

## Print a brief composition summary
pqr

## Examine the storage format (or internal *str*ucture)
str(pqr)

## Print data for the first four atom
pqr$atom[1:4,]

## Print some coordinate data
head(pqr$atom[, c("x","y","z")])

## Print C-alpha coordinates (can also use 'atom.select' function)
head(pqr$atom[pqr$calpha, c("resid","elety","x","y","z")])
inds <- atom.select(pqr, elety="CA")
head( pqr$atom[inds$atom, ] )

## The atom.select() function returns 'indices' (row numbers)
## that can be used for accessing subsets of PDB objects, e.g.
inds <- atom.select(pqr,"ligand")
pqr$atom[inds$atom,]
pqr$xyz[inds$xyz]

## See the help page for atom.select() function for more details.

}
}
\keyword{ IO }
