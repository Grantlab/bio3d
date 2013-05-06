\name{read.pdb}
\alias{read.pdb}
\alias{print.pdb}
\title{ Read PDB File }
\description{
  Read a Protein Data Bank (PDB) coordinate file.
}
\usage{
read.pdb(file, maxlines = 50000, multi = FALSE, rm.insert = FALSE,
         rm.alt = TRUE, het2atom=FALSE, verbose = TRUE)
print.pdb(x, \dots)
}
\arguments{
  \item{file}{ a single element character vector containing the name of
    the PDB file to be read, or the four letter PDB identifier for
    online file access. }
  \item{maxlines}{ the maximum number of lines to read before giving up
    with large files. Default is 50,000 lines. }
  \item{multi}{ logical, if TRUE multiple ATOM records are read for all
    models in multi-model files. }
  \item{rm.insert}{ logical, if TRUE PDB insert records are ignored. }
  \item{rm.alt}{ logical, if TRUE PDB alternate records are ignored. }
  \item{het2atom}{ logical, if TRUE HETATM PDB records are stored as
    ATOM records and returned in the output as such, this should be used
    with caution. }  
  \item{verbose}{ print details of the reading process. }
  \item{x}{ a PDB structure object obtained from
    \code{\link{read.pdb}}.  }
  \item{...}{ additional arguments to \sQuote{print}. }
   
}
\details{
  \code{maxlines} may require increasing for some large multi-model files.
  The preferred means of reading such data is via binary DCD format 
  trajectory files (see the \code{\link{read.dcd}} function).
}
\value{
  Returns a list of class \code{"pdb"} with the following components:
  \item{atom}{ a character matrix containing all atomic coordinate ATOM data, 
  with a row per ATOM and a column per record type.  See below for details
  of the record type naming convention (useful for accessing columns). }
  \item{het }{ a character matrix containing atomic coordinate records
    for atoms within \dQuote{non-standard} HET groups (see \code{atom}). }
  \item{helix }{ \sQuote{start}, \sQuote{end} and \sQuote{length} of H
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{sheet }{ \sQuote{start}, \sQuote{end} and \sQuote{length} of E
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{seqres }{ sequence from SEQRES field. }
  \item{xyz }{ a numeric vector of ATOM coordinate data. }
  \item{xyz.models }{ a numeric matrix of ATOM coordinate data for
    multi-model PDB files. }
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
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
\seealso{ \code{\link{atom.select}}, \code{\link{write.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.fasta.pdb}},
  \code{\link{read.fasta}} }
\examples{
## Read a PDB file from the RCSB online database
pdb <- read.pdb("1bg2")

## Read a PDB file from those included with the package
###pdb <- read.pdb( system.file("examples/1bg2.pdb", package="bio3d") )
  
## Print data for the first atom
pdb$atom[1,]

## Look at the first het atom
pdb$het[1,]

## Print some coordinate data
head(pdb$atom[, c("x","y","z")])

## Or coordinates as a numeric vector
head(pdb$xyz)

## Print C-alpha coordinates (can also use 'atom.select' function)
head(pdb$atom[pdb$calpha, c("resid","elety","x","y","z")])

\dontrun{
## Print SSE data (for helix and sheet),
##  see also dssp and stride functions
pdb$helix
pdb$sheet$start
  
## Print SEQRES data
pdb$seqres

## SEQRES as one leter code
aa321(pdb$seqres)

## Where is the P-loop motif in the ATOM sequence
inds.seq <- motif.find("G....GKT", seq.pdb(pdb))
seq.pdb(pdb)[inds.seq]

## Where is it in the structure
inds.pdb <- atom.select(pdb,resno=inds.seq, elety="CA")
pdb$atom[inds.pdb$atom,]
pdb$xyz[inds.pdb$xyz]


## Renumber residues 
nums <- as.numeric(pdb$atom[,"resno"])
pdb$atom[,"resno"] <- nums - (nums[1] - 1)

## Write out renumbered PDB file
###write.pdb(pdb=pdb,file="eg.pdb")
}
}
\keyword{ IO }
