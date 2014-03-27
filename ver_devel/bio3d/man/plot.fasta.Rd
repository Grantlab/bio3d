\name{plot.fasta}
\alias{plot.fasta}
\title{ Plot a Multiple Sequence Alignment }
\description{
  Produces a schematic representation of a multiple sequence alignment.
}
\usage{
\method{plot}{fasta}(x, plot.labels = TRUE,
                     plot.bars = TRUE,
                     plot.lines = FALSE,
                     plot.axis = TRUE,
                     seq.index = NULL,
                     color.conserved = FALSE,
                     cutoff=0.5, 
                     col=NULL,
                     bars.scale=2,
                     row.spacing=0.5,
                     aln.col="grey50",
                     cex.text=1, add=FALSE, ...)
}
\arguments{
  \item{x}{ multiple sequence alignement of class \sQuote{fasta} as
    obtained from  \code{\link{seqaln}}. }
  \item{plot.labels}{ logical, if TRUE labels will be printed next to
    the sequence bar. }
  \item{plot.bars}{ logical, if TRUE an additional bar representing
    sequence conservation will be plotted. }
  \item{plot.lines}{ logical, if TRUE sequence conservation will be
    represented with a plot. }
  \item{plot.axis}{ logical, if TRUE x-axis will be plotted. }
  \item{seq.index}{ printed tick labels will correspond to the sequence
    of the provided index. }
  \item{color.conserved}{ logical, if TRUE conserved residues will be
    colored according to \dQuote{clustal} coloring scheme. }
  \item{cutoff}{ conservation \sQuote{cutoff} value below which
    alignment columns are not colored. }
  \item{col}{ character vector with color codes for the
    conservation bars. By default, \code{heat.colors} will be used. }
  \item{bars.scale }{ scaling factor for the height of the
    conservation bar when \sQuote{plot.bars=TRUE}. }
  \item{row.spacing }{ space between the sequence bars. }
  \item{aln.col }{ color of the alignment bars. }
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
  \code{\link{entropy}}, \code{\link{aln2html}}. }
\examples{
# Read alignment
aln<-read.fasta(system.file("examples/kif1a.fa",package="bio3d"))

## alignment plot
plot(aln)

\dontrun{
infile <- "http://pfam.sanger.ac.uk/family/PF00071/alignment/seed/format?format=fasta"
aln <- read.fasta( infile )
plot(aln)
}
}
\keyword{ hplot }
