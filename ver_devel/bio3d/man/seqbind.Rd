\name{seqbind}
\alias{seqbind}
\title{ Combine Sequences by Rows Without Recycling }
\description{
  Take vectors and/or matrices arguments and combine them row-wise without
  recycling them (as is the case with \code{\link{rbind}}).
}
\usage{
seqbind(..., blank = "-")
}
\arguments{
  \item{\dots}{ vectors and/or matrices to combine. }
  \item{blank}{ a character to add to short arguments, to achieve the
    same length as the longer argument. }
}
\value{
  A matrix combining input arguments row-wise.
}
\references{
    Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{seqaln}}, \code{\link{read.fasta}},
  \code{\link{read.pdb}}, \code{\link{write.fasta}}, \code{\link{rbind}} }
\examples{
\dontrun{
## Read two pdbs
a.pdb <- read.pdb("1bg2")
b.pdb <- read.pdb("1goj")

seqs <- seqbind(aa321(a.pdb$atom[a.pdb$calpha,"resid"]),
                aa321(b.pdb$atom[b.pdb$calpha,"resid"]))

# seqaln(seqs)
}
}
\keyword{ utilities }
