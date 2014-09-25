\name{write.pdb}
\alias{write.pdb}
\title{ Write PDB Format Coordinate File }
\description{
  Write a Protein Data Bank (PDB) file for a given \sQuote{xyz}
  Cartesian coordinate vector or matrix.
}
\usage{
write.pdb(pdb = NULL, file = "R.pdb", xyz = pdb$xyz, type = NULL, resno = NULL,
 resid = NULL, eleno = NULL, elety = NULL, chain = NULL, insert = NULL, 
 alt = NULL, o = NULL, b = NULL, segid = NULL, elesy = NULL, charge = NULL, 
append = FALSE, verbose = FALSE, chainter = FALSE, end = TRUE, print.segid = FALSE)
}
\arguments{
  \item{pdb}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{file}{ the output file name. }
  \item{xyz}{ Cartesian coordinates as a vector or 3xN matrix. }
  \item{type}{ vector of record types, i.e. "ATOM" or "HETATM", 
    with length equal to length(xyz)/3. }
  \item{resno}{ vector of residue numbers of length equal to
    length(xyz)/3. }
  \item{resid}{ vector of residue types/ids of length equal to
    length(xyz)/3. }
  \item{eleno}{ vector of element/atom numbers of length equal to
    length(xyz)/3. }
  \item{elety}{ vector of element/atom types of length equal to
    length(xyz)/3. }
  \item{chain}{ vector of chain identifiers with length equal to
    length(xyz)/3. }
  \item{insert}{ vector of insertion code with length equal to
    length(xyz)/3. }
  \item{alt}{ vector of alternate record with length equal to
    length(xyz)/3. }
  \item{o}{ vector of occupancy values of length equal to
    length(xyz)/3. }
  \item{b}{ vector of B-factors of length equal to length(xyz)/3. }
  \item{segid}{ vector of segment id of length equal to length(xyz)/3. }
  \item{elesy}{ vector of element symbol of length equal to length(xyz)/3. }
  \item{charge}{ vector of atomic charge of length equal to length(xyz)/3. }
  \item{append}{ logical, if TRUE output is appended to the bottom of an
    existing file (used primarly for writing multi-model files). }
  \item{verbose}{ logical, if TRUE progress details are printed. }
  \item{chainter}{ logical, if TRUE a TER line is inserted at termination
    of a chain. }
  \item{end}{ logical, if TRUE END line is written. }
  \item{print.segid}{ logical, if FALSE segid will not be written. }
}
\details{
  Only the \code{xyz} argument is strictly required.  Other arguments
  assume a default poly-ALA C-alpha structure with a blank chain id,
  occupancy values of 1.00 and B-factors equal to 0.00.
  
  If the input argument \code{xyz} is a matrix then each row is assumed
  to be a different structure/frame to be written to a
  \dQuote{multimodel} PDB file, with frames separated by \dQuote{END}
  records.
  
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of PDB format (version3.3) see:\cr
  \url{http://www.wwpdb.org/documentation/format33/v3.3.html}.
}
\author{
  Barry Grant
  with contributions from Joao Martins.

}
\note{
  Check that:
  (1) \code{chain} is one character long e.g. \dQuote{A}, and 
  (2) \code{resno} and \code{eleno} do not exceed \dQuote{9999}.
}
\seealso{ \code{\link{read.pdb}}, \code{\link{read.dcd}},
  \code{\link{read.fasta.pdb}}, \code{\link{read.fasta}} }
\examples{

# Read a PDB file
pdb <- read.pdb( "1bg2" )

# Renumber residues
nums <- as.numeric(pdb$atom[,"resno"])
nums <- nums - (nums[1] - 1)

# Write out renumbered PDB file
write.pdb(pdb=pdb, resno = nums, file="eg.pdb")

}
\keyword{ IO }
