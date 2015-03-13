\name{as.fasta}
\alias{as.fasta}
\title{ Alignment to FASTA object }
\description{
  Convert alignment/sequence in matrix/vector format to FASTA object.
}
\usage{
as.fasta(x, id=NULL, \dots)
}
\arguments{
  \item{x}{ a sequence character matrix/vector (e.g obtained from
    \code{\link{get.seq}} or \code{\link{seqbind}}). }
  \item{id}{ a vector of sequence names to serve as sequence
    identifers. By default the function will use the row names of the
    alignment if they exists, otherwise ids will be generated. }
  \item{\dots}{ arguments passed to and from functions. }
}
\details{
  This function provides basic functionality to convert a sequence
  character matrix/vector to a FASTA object. 
}
\value{
  Returns a list of class \code{"fasta"} with the following components:
  \item{ali}{ an alignment character matrix with a row per sequence and
    a column per equivalent aminoacid/nucleotide. }
  \item{id}{ sequence names as identifers.}
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{get.seq}}, \code{\link{seqaln}},
  \code{\link{seqbind}}, \code{\link{pdbaln}}
}
\examples{
as.fasta(c("A", "C", "D"))
}
\keyword{ utilities }
