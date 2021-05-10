\name{read.fasta}
\alias{read.fasta}
\title{ Read FASTA formated Sequences }
\description{
  Read aligned or un-aligned sequences from a FASTA format file.
}
\usage{
read.fasta(file, rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
}
\arguments{
  \item{file}{ input sequence file. }
  \item{rm.dup}{ logical, if TRUE duplicate sequences (with the same
    names/ids) will be removed. }
  \item{to.upper}{ logical, if TRUE residues are forced to uppercase. }
  \item{to.dash}{ logical, if TRUE \sQuote{.} gap characters are
    converted to \sQuote{-} gap characters. }  
}
\value{
  A list with two components:
  \item{ali }{ an alignment character matrix with a row per sequence and
    a column per equivalent aminoacid/nucleotide. }
  \item{ids }{ sequence names as identifers.}
  \item{call}{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{  
  For a description of FASTA format see:
  \url{https://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml}.
  When reading alignment files, the dash \sQuote{-} is interpreted as
  the gap character.
}
\seealso{ \code{\link{read.fasta.pdb}} }
\examples{
# Read alignment
aln <- read.fasta(system.file("examples/hivp_xray.fa",package="bio3d"))

# Print alignment overview
aln

# Sequence names/ids
head( aln$id )

# Alignment positions 335 to 339
head( aln$ali[,33:39] )

# Sequence d2a4f_b
aa123( aln$ali["d2a4f_b",] )

# Write out positions 33 to 45 only
#aln$ali=aln$ali[,30:45]
#write.fasta(aln, file="eg2.fa")

}
\keyword{ IO }
