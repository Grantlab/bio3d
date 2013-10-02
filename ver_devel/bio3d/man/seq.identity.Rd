\name{seq.identity}
\alias{seq.identity}
\title{ Percent Identity}
\description{
  Determine the percent identity scores for aligned sequences.
}
\usage{
  seq.identity(alignment, normalize=TRUE, ncore=1, nseg.scale=1)
}
\arguments{
  \item{alignment}{ sequence alignment obtained from
    \code{\link{read.fasta}} or an alignment character matrix. }
  \item{normalize}{ logical, if TRUE output is normalized to values
    between 0 and 1 otherwise percent identity is returned. }
  \item{ncore}{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. } 
  \item{nseg.scale}{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
  The percent identity value is a single numeric score determined for
  each pair of aligned sequences.  It measures the number of identical
  residues (\dQuote{matches}) in relation to the length of the alignment.
}
\value{
  Returns a numeric matrix with all pairwise identity values.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.fasta}}, \code{\link{ide.filter}}, \code{\link{entropy}}, \code{\link{consensus}} }
\examples{
\dontrun{
#aln <- read.fasta( system.file("examples/kinesin_xray.fa",
#                     package = "bio3d") )
}

data(kinesin)
attach(kinesin, warn.conflicts=FALSE)

ide.mat <- seq.identity(aln)

# Plot identity matrix
plot.dmat(ide.mat, color.palette=mono.colors,
          main="Sequence Identity", xlab="Structure No.",
          ylab="Structure No.")

# Histogram of pairwise identity values
hist(ide.mat[upper.tri(ide.mat)], breaks=30,xlim=c(0,1),
     main="Sequence Identity", xlab="Identity")

# Compare two sequences
seq.identity( rbind(aln$ali[1,], aln$ali[20,]) )

}
\keyword{ utilities }
