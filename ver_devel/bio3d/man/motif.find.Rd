\name{motif.find}
\alias{motif.find}
\title{ Find Sequence Motifs. }
\description{
  Return Position Indices of a Short Sequence Motif Within a Larger Sequence.
}
\usage{
motif.find(motif, sequence)
}
\arguments{
  \item{motif}{ a character vector of the short sequence motif. }
  \item{sequence}{ a character vector of the larger sequence. }
}
\details{
  The sequence and the motif can be given as a either a multiple or
  single element character vector. The dot character and other valid
  \code{regexpr} characters are allowed in the motif, see examples.
}
\value{
  Returns a vector of position indices within the sequence where the
  motif was found, see examples.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{regexpr}}, \code{\link{read.fasta}}, \code{\link{pdbseq}} 
}
\examples{
\donttest{
aa.seq <- pdbseq( read.pdb( get.pdb("4q21", URLonly=TRUE) ) )
motif = c("G....GKS")
motif.find(motif, aa.seq)
}
}
\keyword{ utilities }
