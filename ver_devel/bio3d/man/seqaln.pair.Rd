\name{seqaln.pair}
\alias{seqaln.pair}
\title{ Sequence Alignment of Identical Protein Sequences }
\description{
  Create multiple alignments of amino acid sequences according to the
  method of Edgar.
}
\usage{
seqaln.pair(aln, extra.args = "", ...)
}
\arguments{
  \item{aln}{ a sequence character matrix, as obtained from
    \code{\link{seqbind}}, or an alignment list object as obtained from
      \code{\link{read.fasta}}. }
  \item{extra.args}{ a single character string containing extra command
    line arguments for the alignment program. }
  \item{\dots}{ additional arguments for the function \code{\link{seqaln}}. }
}
\details{
  This function is intended for the alignment of identical sequences only.
  For standard alignment see the related function \code{\link{seqaln}}.

  This function is useful for determining the equivalences between
  sequences and structures. For example in aligning a PDB sequence to an
  existing multiple sequence alignment, where one would first mask the
  alignment sequences and then run the alignment to determine equivalences.
}
\value{
  A list with two components:
  \item{ali}{ an alignment character matrix with a row per sequence and
    a column per equivalent aminoacid/nucleotide. }
  \item{ids}{ sequence names as identifers.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \sQuote{MUSCLE} is the work of Edgar:
  Edgar (2004) \emph{Nuc. Acid. Res.} \bold{32}, 1792--1797.

  Full details of the \sQuote{MUSCLE} algorithm, along with download and
  installation instructions can be obtained from:\cr
  \url{http://www.drive5.com/muscle}.
}
\author{ Barry Grant }
\note{
  A system call is made to the \sQuote{MUSCLE} program, which must be
  installed on your system and in the search path for executables.
}
\seealso{ \code{\link{seqaln}}, \code{\link{read.fasta}},
  \code{\link{read.fasta.pdb}}, \code{\link{seqbind}}
}
\examples{

## NOTE: FOLLOWING EXAMPLE NEEDS MUSCLE INSTALLED
if(check.utility("muscle")) {

##- Aligning a PDB sequence to an existing sequence alignment


##- Simple example
aln <- list(ali=seqbind(c("X","C","X","X","A","G","K"),
                        c("C","-","A","X","G","X","X","K")),
            id=c("a","b"))

seqaln.pair(aln)

}
}
\keyword{ utilities }
