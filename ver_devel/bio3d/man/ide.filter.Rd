\name{ide.filter}
\alias{ide.filter}
\title{ Percent Identity Filter }
\description{
  Identify and filter subsets of sequences at a given sequence identity cutoff.
}
\usage{
ide.filter(aln = NULL, ide = NULL, cutoff = 0.6, verbose = TRUE, ncore=1, nseg.scale=1)
}
\arguments{
  \item{aln}{ sequence alignment list, obtained from
    \code{\link{seqaln}} or \code{\link{read.fasta}}, or an alignment
    character matrix. Not used if \sQuote{ide} is given.}
  \item{ide}{ an optional identity matrix obtained from
    \code{\link{seqidentity}}. }
  \item{cutoff}{ a numeric identity cutoff value ranging between 0 and 1. }
  \item{verbose}{ logical, if TRUE print details of the clustering process. }
  \item{ncore}{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. } 
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
  This function performs hierarchical cluster analysis of a given sequence
  identity matrix \sQuote{ide}, or the identity matrix calculated from a given
  alignment \sQuote{aln}, to identify sequences that fall below a given
  identity cutoff value \sQuote{cutoff}.  
}
\value{
  Returns a list object with components:
  \item{ind}{indices of the sequences below the cutoff value.}
  \item{tree}{an object of class \code{"hclust"}, which describes the
    tree produced by the clustering process. }
  \item{ide}{a numeric matrix with all pairwise identity values.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.fasta}}, \code{\link{seqaln}},
  \code{\link{seqidentity}}, \code{\link{entropy}}, \code{\link{consensus}} }
\examples{
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)
ide.mat <- seqidentity(aln)

# Histogram of pairwise identity values
par(mfrow=c(2,1))
hist(ide.mat[upper.tri(ide.mat)], breaks=30,xlim=c(0,1),
     main="Sequence Identity", xlab="Identity")

k <- ide.filter(ide=ide.mat, cutoff=0.6)
ide.cut <- seqidentity(aln$ali[k$ind,])
hist(ide.cut[upper.tri(ide.cut)], breaks=10, xlim=c(0,1),
     main="Sequence Identity", xlab="Identity")

#plot(k$tree, axes = FALSE, ylab="Sequence Identity")
#print(k$ind) # selected
}
\keyword{ utilities }
