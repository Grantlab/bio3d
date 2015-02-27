\name{read.pdb}
\alias{read.pdb}
\alias{print.pdb}
\alias{summary.pdb}
\title{ Read PDB File}
\description{
  Read a Protein Data Bank (PDB) coordinate file.
}
\usage{
read.pdb(file, maxlines = -1, multi = FALSE, rm.insert = FALSE,
         rm.alt = TRUE, ATOM.only = FALSE, verbose = TRUE)

\method{print}{pdb}(x, printseq=TRUE, \dots)

\method{summary}{pdb}(object, printseq=FALSE, \dots)
}
\arguments{
  \item{file}{ a single element character vector containing the name of
    the PDB file to be read, or the four letter PDB identifier for
    online file access. }
  \item{maxlines}{ the maximum number of lines to read before giving up
    with large files. By default if will read up to the end of input on
    the connection. }
  \item{multi}{ logical, if TRUE multiple ATOM records are read for all
    models in multi-model files and their coordinates returned. }
  \item{rm.insert}{ logical, if TRUE PDB insert records are ignored. }
  \item{rm.alt}{ logical, if TRUE PDB alternate records are ignored. }
  \item{ATOM.only}{ logical, if TRUE only ATOM/HETATM records are stored.
    Useful for speed enhancements with large files where secondary structure, 
    biological unit and other remark records are not required. }
  \item{verbose}{ print details of the reading process. }
  \item{x}{ a PDB structure object obtained from
    \code{\link{read.pdb}}.  }
  \item{object}{ a PDB structure object obtained from
    \code{\link{read.pdb}}.  }
  \item{printseq}{ logical, if TRUE the PDB ATOM sequence will be printed
    to the screen. See also \code{\link{pdbseq}}.  }
  \item{...}{ additional arguments to \sQuote{print}. }
   
}
\details{
  \code{maxlines} may be set so as to restrict the reading to a portion
  of input files. Note that the preferred means of reading large 
  multi-model files is via binary DCD or NetCDF format trajectory files 
  (see the \code{\link{read.dcd}} and \code{\link{read.ncdf}} functions).
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
  \item{remark }{ a list object containing information taken from 'REMARK'
    records of a \code{"pdb"}. It can be used for building biological units
    (See \code{\link{biounit}}). }
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
\seealso{
  \code{\link{atom.select}}, \code{\link{write.pdb}},
  \code{\link{trim.pdb}}, \code{\link{cat.pdb}}, 
  \code{\link{read.prmtop}}, \code{\link{as.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.ncdf}}, 
  \code{\link{read.fasta.pdb}}, \code{\link{read.fasta}}, 
  \code{\link{biounit}}
}
\examples{
## Read a PDB file from the RCSB online database
#pdb <- read.pdb("4q21")

## Read a PDB file from those included with the package
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Print a brief composition summary
pdb

## Examine the storage format (or internal *str*ucture)
str(pdb)

## Print data for the first four atom
pdb$atom[1:4,]

## Print some coordinate data
head(pdb$atom[, c("x","y","z")])

## Or coordinates as a numeric vector
head(pdb$xyz)

## Print C-alpha coordinates (can also use 'atom.select' function)
head(pdb$atom[pdb$calpha, c("resid","elety","x","y","z")])
inds <- atom.select(pdb, elety="CA")
head( pdb$atom[inds$atom, ] )

## The atom.select() function returns 'indices' (row numbers)
## that can be used for accessing subsets of PDB objects, e.g.
inds <- atom.select(pdb,"ligand")
pdb$atom[inds$atom,]
pdb$xyz[inds$xyz]

## See the help page for atom.select() function for more details.


\dontrun{
## Print SSE data for helix and sheet,
##  see also dssp() and stride() functions
print.sse(pdb)
pdb$helix
pdb$sheet$start
  
## Print SEQRES data
pdb$seqres

## SEQRES as one letter code
aa321(pdb$seqres)

## Where is the P-loop motif in the ATOM sequence
inds.seq <- motif.find("G....GKT", pdbseq(pdb))
pdbseq(pdb)[inds.seq]

## Where is it in the structure
inds.pdb <- atom.select(pdb,resno=inds.seq, elety="CA")
pdb$atom[inds.pdb$atom,]
pdb$xyz[inds.pdb$xyz]

## View in interactive 3D mode
#view(pdb)
}

}
\keyword{ IO }

