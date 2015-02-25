\name{read.crd.charmm}
\alias{read.crd.charmm}
\title{ Read CRD File }
\description{
  Read a CHARMM CARD (CRD) coordinate file.
}
\usage{
\method{read.crd}{charmm}(file, ext = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{file}{ the name of the CRD file to be read.  }
  \item{ext}{logical, if TRUE assume expanded CRD format. }
  \item{verbose}{ print details of the reading process. }
  \item{\dots}{ arguments going nowhere. }
}
\details{
  See the function \code{\link{read.pdb}} for more details.
}
\value{
  Returns a list with the following components:
  \item{atom}{ a character matrix containing all atomic coordinate data,
  with a row per atom and a column per record type.  See below for details
  of the record type naming convention (useful for accessing columns). }
  \item{xyz }{ a numeric vector of coordinate data. }
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  For a description of CHARMM CARD (CRD) format see:\cr
  \url{http://www.charmmtutorial.org/index.php/CHARMM:The_Basics}.
}
\author{ Barry Grant }
\note{
  Similar to the output of \code{\link{read.pdb}}, the column names of
  \code{atom} can be used as a convenient means of data access, namely:
  Atom serial number \dQuote{eleno},
  Atom type \dQuote{elety},
  Alternate location indicator \dQuote{alt},
  Residue name \dQuote{resid},
  Residue sequence number \dQuote{resno},
  Code for insertion of residues \dQuote{insert},
  Orthogonal coordinates \dQuote{x},
  Orthogonal coordinates \dQuote{y},
  Orthogonal coordinates \dQuote{z},
  Weighting factor \dQuote{b}.
  See examples for further details.
}
\seealso{ \code{\link{write.crd}}, \code{\link{read.pdb}},
  \code{\link{atom.select}}, \code{\link{write.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.fasta.pdb}},
  \code{\link{read.fasta}} }
\examples{
\dontrun{
 pdb <- read.pdb("1bg2")
 crdfile <- tempfile()
 write.crd(pdb, file=crdfile)
 crd <- read.crd(crdfile)
 ca.inds <- which(crd$calpha)
 crd$atom[ca.inds[1:20],c("x","y","z")]
# write.pdb(crd, file=tempfile())
}
}
\keyword{ IO }
