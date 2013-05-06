\name{write.pqr}
\alias{write.pqr}
\title{ Write PQR Format Coordinate File }
\description{
  Write a PQR file for a given \sQuote{xyz} Cartesian coordinate vector
  or matrix.
}
\usage{
write.pqr(pdb = NULL, xyz = pdb$xyz, resno = NULL, resid = NULL, eleno =
NULL, elety = NULL, chain = NULL, o = NULL, b = NULL, het = FALSE,
append = FALSE, verbose = FALSE, chainter = FALSE, file = "R.pdb")
}
\arguments{
  \item{pdb}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{xyz}{ Cartesian coordinates as a vector or 3xN matrix. }
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
  \item{o}{ vector of occupancy values of length equal to
    length(xyz)/3. }
  \item{b}{ vector of B-factors of length equal to length(xyz)/3. }
  \item{het}{ logical, if TRUE \sQuote{HETATM} records from \code{pdb}
    object are written to the output PDB file. } 
  \item{append}{ logical, if TRUE output is appended to the bottom of an
    existing file (used primarly for writing multi-model files). }
  \item{verbose}{ logical, if TRUE progress details are printed. }
  \item{chainter}{ logical, if TRUE a TER line is inserted between
    chains. }
  \item{file}{ the output file name. }
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
pdb <- read.pdb(system.file("examples/d1bg2__.ent", package="bio3d"))

# Renumber residues
nums <- as.numeric(pdb$atom[,"resno"])
nums <- nums - (nums[1] - 1)

# Write out renumbered PDB file
write.pdb(pdb=pdb, resno = nums, file="eg.pdb")

}
\keyword{ IO }
