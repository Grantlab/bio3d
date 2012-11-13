\name{gap.inspect}
\alias{gap.inspect}
\title{ Alignment Gap Summary }
\description{
  Report the number of gaps per sequence and per position for a given alignment.
}
\usage{
gap.inspect(x)
}
\arguments{
  \item{x}{ a matrix or an alignment data structure obtained from
    \code{\link{read.fasta}} or \code{\link{read.fasta.pdb}}.}
}
\details{
  Reports the number of gap characters per row (i.e. sequence) and
  per column (i.e. position) for a given \code{alignment}. In addition,
  the indices for gap and non-gap containing coloums are returned along
  with a binary matrix indicating the location of gap positions. 
}
\value{
  Returns a list object with the following components:
  \item{row }{a numeric vector detailing the number of gaps per row
    (i.e. sequence).}
  \item{col }{a numeric vector detailing the number of gaps per column
    (i.e. position).}
  \item{t.inds }{ indices for gap containing coloums }
  \item{f.inds }{ indices for non-gap containing coloums}  
  \item{bin }{a binary numeric matrix with the same dimensions as the
    \code{alignment}, with 0 at non-gap positions and 1 at gap
    positions.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  During alignment, gaps are introduced into sequences that are believed
  to have undergone deletions or insertions with respect to other
  sequences in the alignment.  These gaps, often referred to as
  indels, can be represented with \sQuote{NA}, a \sQuote{-} or
  \sQuote{.} character.

  This function gives an overview of gap occurrence and may be useful
  when considering positions or sequences that could/should be excluded
  from further analysis.
}
\seealso{ \code{\link{read.fasta}}, \code{\link{read.fasta.pdb}}}
\examples{
aln <- read.fasta( system.file("examples/hivp_xray.fa",
                   package = "bio3d") )

gap.stats <- gap.inspect(aln$ali)
gap.stats$row # Gaps per sequence
gap.stats$col # Gaps per position
##gap.stats$bin # Binary matrix (1 for gap, 0 for aminoacid)
##aln[,gap.stats$f.inds] # Alignment without gap positions 

plot(gap.stats$col, typ="h", ylab="No. of Gaps") 
}
\keyword{ utilities }
