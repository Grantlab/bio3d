\name{plot.fasta}
\alias{plot.fasta}
\title{ Plot a Multiple Sequence Alignment }
\description{
  Produces a schematic representation of a multiple sequence alignment.
}
\usage{
\method{plot}{fasta}(x, hc = TRUE, labels = x$id, cex.lab = 0.7,
                         xlab = "Alignment index",
                         main = "Sequence Alignment Overview",
                         mar4 = 4, ...)
}
\arguments{
  \item{x}{ multiple sequence alignement of class \sQuote{fasta} as
    obtained from  \code{\link{seqaln}}. }
  \item{hc}{ logical, if TRUE plot a dendrogram on the left
    side. Alternatively, an object obtained from \code{\link{hclust}}
    can be provided. }
  \item{labels}{ labels corresponding to each row in the alignment. }
  \item{cex.lab}{ scaling factor for the labels. }
  \item{xlab}{ label for x-axis. }
  \item{main}{ a main title for the plot. }
  \item{mar4}{ margin size for the labels. }
  \item{\dots}{ additional arguments passed to function \code{\link{hclust}}. }
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
  \code{\link{entropy}}, \code{\link{aln2html}}. }
\examples{
# Read alignment
aln <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))

## alignment plot
plot(aln, labels=basename.pdb(aln$id))

## Works also for a 'pdbs' object
attach(transducin)
plot(pdbs)

\dontrun{
infile <- "http://pfam.xfam.org/family/PF00071/alignment/seed/format?format=fasta"
aln <- read.fasta( infile )
plot(aln)
}
}
\keyword{ hplot }
