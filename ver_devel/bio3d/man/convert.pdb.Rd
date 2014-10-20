\name{convert.pdb}
\alias{convert.pdb}
\title{Renumber and Convert Between Various PDB formats}
\description{
  Renumber and convert between CHARMM, Amber, Gromacs and Brookhaven
  PDB formats.
}
\usage{
convert.pdb(pdb, type="original",   
                renumber = FALSE, first.resno = 1, first.eleno = 1,
                consecutive=TRUE, rm.h = TRUE, rm.wat = FALSE,
                verbose=TRUE)
}
\arguments{
  \item{pdb}{ a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{type}{ output format, one of \sQuote{original}, \sQuote{pdb}, 
    \sQuote{charmm}, \sQuote{amber}, or \sQuote{gromacs}. The default 
    option of \sQuote{original} results in no conversion. }
  \item{renumber}{ logical, if TRUE atom and residue records are
    renumbered using \sQuote{first.resno} and \sQuote{first.eleno}. }
  \item{first.resno}{ first residue number to be used if
    \sQuote{renumber} is TRUE. }
  \item{first.eleno}{ first element number to be used if
    \sQuote{renumber} is TRUE. }
  \item{consecutive}{ logical, if TRUE renumbering will result in 
    consecutive residue numbers spanning all chains. Otherwise new residue  
    numbers will begin at \sQuote{first.resno} for each chain. }
  \item{rm.h}{ logical, if TRUE hydrogen atoms are removed. }
  \item{rm.wat}{ logical, if TRUE water atoms are removed. }
  \item{verbose}{ logical, if TRUE details of the conversion process 
    are printed. }
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
\donttest{

# Read a PDB file
pdb <- read.pdb("4q21")
pdb
head( pdb$atom[pdb$calpha,"resno"] )

# Convert to CHARMM format
new <- convert.pdb(pdb, type="amber", renumber=TRUE, first.resno=22 )
head( new$atom[new$calpha,"resno"] )

# Write a PDB file
#write.pdb(new, file="tmp4amber.pdb")

}
}
\keyword{ utilities }
