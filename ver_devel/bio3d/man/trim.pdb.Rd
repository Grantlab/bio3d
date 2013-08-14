\name{trim.pdb}
\alias{trim.pdb}
\title{ Trim a PDB Object To A Subset of Atoms. }
\description{
  Produce a new smaller PDB object, containing a subset of atoms, from a
  given larger PDB object. 
}
\usage{
trim.pdb(pdb, inds = NULL)
}
\arguments{
  \item{pdb}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{inds}{ a list object of ATOM and XYZ indices as obtained from
    \code{\link{atom.select}}. }
}
\details{
  This is a basic utility function for creating a new PDB object based
  on a selection of atoms.
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
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of PDB format (version3.3) see:\cr
  \url{http://www.wwpdb.org/documentation/format33/v3.3.html}.
.
}
\author{ Barry Grant, Lars Skjaerven }
\note{
  \code{het}, \code{helix}, \code{sheet}, and \code{seqres} list
  components are returned unmodified.
  
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
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}},
  \code{\link{write.pdb}} }
\examples{
\dontrun{
## Read a PDB file from the RCSB online database
pdb <- read.pdb("1bg2")

## Select calpha atoms
sele <- atom.select(pdb, "calpha")

## Trim PDB
new.pdb <- trim.pdb(pdb, inds=sele)

## Write to file
write.pdb(new.pdb, file="calpha.pdb")
}
}
\keyword{ utilities }
