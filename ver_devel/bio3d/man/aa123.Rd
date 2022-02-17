\name{aa123}
\alias{aa123}
\alias{aa321}
\title{ Convert Between 1-letter and 3-letter Aminoacid Codes }
\description{
  Convert between one-letter IUPAC aminoacid codes and three-letter PDB
  style aminoacid codes.
}
\usage{
aa123(aa)
aa321(aa)
}
\arguments{
  \item{aa}{ a character vector of individual aminoacid codes. }
}
\details{
  Standard conversions will map \sQuote{A} to \sQuote{ALA}, \sQuote{G} to
  \sQuote{GLY}, etc.
  Non-standard codes in \code{aa} will generate a warning and return
  \sQuote{UNK} or \sQuote{X}.
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
\seealso{ \code{\link{read.pdb}}, \code{\link{read.fasta}}, \code{\link{pdbseq}} }
\examples{
# Simple conversion
aa123(c("D","L","A","G","S","H"))
aa321(c("ASP", "LEU", "ALA", "GLY", "SER", "HIS"))

\dontrun{
# Extract sequence from a PDB file's ATOM and SEQRES cards
pdb <- read.pdb("1BG2") 
s <- aa321(pdb$seqres)                   # SEQRES
a <- aa321(pdb$atom[pdb$calpha,"resid"]) # ATOM

# Write both sequences to a fasta file
write.fasta(alignment=seqbind(s,a), id=c("seqres","atom"), file="eg2.fa")

# Alternative approach for ATOM sequence extraction
pdbseq(pdb)
pdbseq(pdb, aa1=FALSE )
}
}
\keyword{ utilities }
