\name{write.crd}
\alias{write.crd}
\title{ Write CRD File }
\description{
  Write a CHARMM CARD (CRD) coordinate file.
}
\usage{
write.crd(pdb = NULL, xyz = pdb$xyz, resno = NULL, resid = NULL, 
   eleno = NULL, elety = NULL, segid = NULL, resno2 = NULL, b = NULL, 
   verbose = FALSE, file = "R.crd")
}
\arguments{
  \item{pdb}{ a structure object obtained from \code{\link{read.pdb}} or
    \code{\link{read.crd}}. }
  \item{xyz}{ Cartesian coordinates as a vector or 3xN matrix. }
  \item{resno}{ vector of residue numbers of length equal to
    length(xyz)/3. }
  \item{resid}{ vector of residue types/ids of length equal to
    length(xyz)/3. }
  \item{eleno}{ vector of element/atom numbers of length equal to
    length(xyz)/3. }
  \item{elety}{ vector of element/atom types of length equal to
    length(xyz)/3. }
  \item{segid}{ vector of segment identifiers with length equal to
    length(xyz)/3. }
  \item{resno2}{ vector of alternate residue numbers of length equal to
    length(xyz)/3. }
  \item{b}{ vector of weighting factors of length equal to length(xyz)/3. }
  \item{verbose}{ logical, if TRUE progress details are printed. }
  \item{file}{ the output file name. }
}
\details{
  Only the \code{xyz} argument is strictly required.  Other arguments
  assume a default poly-ALA C-alpha structure with a blank segid and
  B-factors equal to 0.00.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  For a description of CHARMM CARD (CRD) format see:\cr
  \url{https://www.charmm.org/wiki//index.php/CHARMM:The_Basics}.
}
\author{ Barry Grant }
\note{
  Check that \code{resno} and \code{eleno} do not exceed \dQuote{9999}.
}
\seealso{ \code{\link{read.crd}}, \code{\link{read.pdb}},
  \code{\link{atom.select}}, \code{\link{write.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.fasta.pdb}},
  \code{\link{read.fasta}}
}
\examples{
\dontrun{
# Read a PDB file
pdb <- read.pdb( "1bg2" )
summary(pdb)
# Convert to CHARMM format
new <- convert.pdb(pdb, type="charmm")
summary(new)
# Write a CRD file
write.crd(new, file="4charmm.crd")
}
}
\keyword{ IO }
