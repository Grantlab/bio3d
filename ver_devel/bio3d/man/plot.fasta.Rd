\name{plot.fasta}
\alias{plot.fasta}
\title{ Plot Multiple Sequence Alignment }
\description{
  Produces a schematic representation of a multiple sequence alignment.
}
\usage{
\method{plot}{fasta}(x, plot.labels=TRUE, plot.conserv=TRUE,
                     conserv.col=NULL, conserv.height.scale=1,
                     conserv.intervals=c(1, .75, .5, .25),
                     row.height=100, row.spacer=20,
                     cex.text=1, add=FALSE, ...)
}
\arguments{
  \item{x}{ multiple sequence alignement of class \sQuote{fasta} as
    obtained from  \code{\link{seqaln}}. }
  \item{plot.labels}{ logical, if TRUE labels will be printed next to
    the sequence bar. }
  \item{plot.conserv}{ logical, if TRUE an additional bar representing
    sequence conservation will be plotted. }
  \item{conserv.col}{ character vector with color codes for the
    conservation bar. Must be of length(conserv.intervals) -
    1. By default, \code{heat.colors} will be used. }
  \item{conserv.height.scale }{ scaling factor for the height of the
    conservation bar. }
  \item{conserv.intervals }{ numeric vector providing the intervals at
    which the conservation bar should be plotted. }
  \item{row.height }{ height of each sequence bar. } 
  \item{row.spacer }{ space between the sequence bars. }
  \item{cex.text }{ scaling factor for the labels. }
  \item{add }{ logical, if TRUE \code{plot.new()} will not be called. }
  \item{\dots}{ additional arguments not paseed anywhere. }
}
\details{
  \code{plot.fasta} is a utility function for producting a schematic
  representation of a multiple sequence alignment.  
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{seqaln}}, \code{\link{read.fasta}},
  \code{\link{entropy}}. }
\examples{
# Read alignment
aln<-read.fasta(system.file("examples/kif1a.fa",package="bio3d"))

## alignment plot
plot(aln)


\dontrun{
seq <- get.seq("P00439")
blast <- blast.pdb(paste(seq$ali, collapse=""), database="swissprot")

hits <- plot.blast(blast, cutoff = 200)
all.seqs <- get.seq(hits$gi.id)
aln <- seqaln(all.seqs)

plot(aln, cex.text=0.5, conserv.height.scale=3)
}
}
\keyword{ hplot }
