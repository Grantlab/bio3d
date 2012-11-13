\name{convert.pdb}
\alias{convert.pdb}
\title{Convert Between Various PDB formats}
\description{
  Convert between CHARMM, Amber, Gromacs and Brookhaven PDB
  formats.
}
\usage{
convert.pdb(pdb, type, renumber = FALSE, first.resno = 1, first.eleno = 1, rm.h = TRUE, rm.wat = FALSE)
}
\arguments{
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{type}{ output format. }
  \item{renumber}{ logical, if TRUE atom and residue records are
    renumbered using \sQuote{first.resno} and \sQuote{first.eleno}. }
  \item{first.resno}{ first residue number to be used if
    \sQuote{renumber} is TRUE. }
  \item{first.eleno}{ first element number to be used if
    \sQuote{renumber} is TRUE. }
  \item{rm.h}{ logical, if TRUE hydrogen atoms are removed. }
  \item{rm.wat}{ logical, if TRUE water atoms are removed. }
}
\details{
  Convert atom names and residue names, renumber atom and residue
  records, strip water and hydrogen atoms from \code{pdb} objects.

  Format \code{type} can be one of \dQuote{ori}, \dQuote{pdb}, \dQuote{charmm},
      \dQuote{amber} or \dQuote{gromacs}.
}
\value{
  Returns a list of class \code{"pdb"}, with the following components:
  \item{atom}{ a character matrix containing all atomic coordinate ATOM data,
  with a row per ATOM and a column per record type.  See below for details
  of the record type naming convention (useful for accessing columns). }
  \item{het}{ a character matrix containing atomic coordinate records
    for atoms within \dQuote{non-standard} HET groups (see \code{atom}). }
  \item{helix}{ \sQuote{start}, \sQuote{end} and \sQuote{length} of H
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{sheet}{ \sQuote{start}, \sQuote{end} and \sQuote{length} of E
    type sse, where start and end are residue numbers \dQuote{resno}. }
  \item{seqres }{ sequence from SEQRES field. }
  \item{xyz}{ a numeric vector of ATOM coordinate data. }
  \item{calpha}{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  For a description of PDB format (version2.2) see:\cr
  \url{http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html}.
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
# Read a PDB file
pdb <- read.pdb( system.file("examples/1bg2.pdb", package="bio3d") )
pdb.summary(pdb)
# Convert to CHARMM format
new <- convert.pdb(pdb, type="charmm")
pdb.summary(new)
# Write a PDB and CRD file
#write.pdb(new, file="4charmm.pdb")
#write.crd(new, file="4charmm.crd")
}
\keyword{ utilities }
