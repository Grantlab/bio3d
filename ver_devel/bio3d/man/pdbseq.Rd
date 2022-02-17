\name{pdbseq}
\alias{pdbseq}
\title{ Extract The Aminoacid Sequence From A PDB Object }
\description{
  Return a vector of the one-letter IUPAC or three-letter PDB
  style aminoacid codes from a given PDB object.
}
\usage{
pdbseq(pdb, inds = NULL, aa1 = TRUE)
}
\arguments{
  \item{pdb}{a PDB structure object obtained from
    \code{\link{read.pdb}}. }
  \item{inds}{ a list object of ATOM and XYZ indices as obtained from
    \code{\link{atom.select}}. }
  \item{aa1}{ logical, if TRUE then the one-letter IUPAC sequence is
    returned. IF FALSE then the three-letter PDB style sequence is
    returned.}
}
\details{
  See the examples below and the functions \code{\link{atom.select}} and \code{\link{aa321}}
  for further details.
}
\value{
   A character vector of aminoacid codes.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
  For a description of IUPAC one-letter codes see:\cr
  \url{https://www.insdc.org/documents/feature_table.html#7.4.3}

  For more information on PDB residue codes see:\cr
  \url{http://ligand-expo.rcsb.org/ld-search.html}

}
\author{ Barry Grant }
\seealso{ \code{\link{read.pdb}}, \code{\link{atom.select}},
  \code{\link{aa321}}, \code{\link{read.fasta}} }
\examples{
\dontrun{
pdb <- read.pdb( "5p21" )
pdbseq(pdb)

#pdbseq(pdb, inds=atom.select(pdb, resno=5:15, elety="CA"), aa1=FALSE)
}
}
\keyword{ utilities }
