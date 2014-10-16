\name{is.gap}
\alias{is.gap}
\title{ Gap Characters }
\description{
  Test for the presence of gap characters.
}
\usage{
is.gap(x, gap.char = c("-", "."))
}
\arguments{
  \item{x}{ an R object to be tested. Typically a sequence vector or 
  sequence/structure alignment object as returned from \code{seqaln}, 
  \code{pdbaln} etc. }
  \item{gap.char}{ a character vector containing the gap character types
    to test for. }
}
\value{
  Returns a logical vector with the same length as the input vector, or 
  the same length as the number of columns present in an alignment input 
  object \sQuote{x}. In the later case TRUE elements corresponding to 
  \sQuote{gap.char} matches in any alignment column (i.e. gap containing 
  columns).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  During alignment, gaps are introduced into sequences that are believed
  to have undergone deletions or insertions with respect to other
  sequences in the alignment.  These gaps, often referred to as
  indels, can be represented with \sQuote{NA}, \sQuote{-} or
  \sQuote{.} characters.

  This function provides a simple test for the presence of such
  characters, or indeed any set of user defined characters set by
  the \sQuote{gap.char} argument.
}
\seealso{ \code{\link{gap.inspect}}, \code{\link{read.fasta}},
  \code{\link{read.fasta.pdb}}, \code{\link{seqaln}}, \code{\link{pdbaln}} }
\examples{

is.gap( c("G",".","X","-","G","K","S","T") )

\dontrun{
aln <- read.fasta( system.file("examples/kif1a.fa",
                   package = "bio3d") )

##- Print only non-gap positions (i.e. no gaps in any sequence)
aln$ali[, !is.gap(aln) ]

##- Mask any existing gaps with an "X"
xaln <- aln
xaln$ali[ is.gap(xaln$ali) ]="X"

##- Read a new PDB and align its sequence to the existing masked alignment
pdb <- read.pdb( "1mkj" )
seq2aln(pdbseq(pdb), xaln, id = "1mkj")
}
}
\keyword{ utilities }
