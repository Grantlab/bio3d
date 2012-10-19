\name{consensus}
\alias{consensus}
\title{ Sequence Consensus for an Alignment }
\description{
  Determines the consensus sequence for a given alignment at a given
  identity cutoff value.
}
\usage{
  consensus(alignment, cutoff = 0.6)
}
\arguments{
  \item{alignment}{ an \code{alignment} object created by the
    \code{\link{read.fasta}} function or an alignment character matrix. }
  \item{cutoff}{ a numeric value beteen 0 and 1, indicating the minimum
    sequence identity threshold for determining a consensus amino acid.
    Default is 0.6, or 60 percent residue identity. }
}
\value{
  A vector containing the consensus sequence, where \sQuote{-} represents
  positions with no consensus (i.e. under the \code{cutoff})
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{read.fasta}} }
\examples{

#-- Read kinesin alignment
aln <- read.fasta(system.file("examples/kinesin_xray.fa",package="bio3d"))

# Generate consensus
con <- consensus(aln)
print(con$seq)

# Generate consensus for sub-sequence
con <- consensus(aln$ali[,330:474])
print(con$seq)


#-- Read HIV protease alignment
aln <- read.fasta(system.file("examples/hivp_xray.fa",package="bio3d"))
con <- consensus(aln$ali)

# Plot residue frequency matrix
##png(filename = "freq.png", width = 1500, height = 780)
col <- mono.colors(32)
aa  <- rev(rownames(con$freq))

image(x=1:ncol(con$freq),
      y=1:nrow(con$freq),
      z=as.matrix(rev(as.data.frame(t(con$freq)))),
      col=col, yaxt="n", xaxt="n",
      xlab="Alignment Position", ylab="Residue Type")

# Add consensus along the axis
axis(side=1, at=seq(0,length(con$seq),by=5))
axis(side=2, at=c(1:22), labels=aa)
axis(side=3, at=c(1:length(con$seq)), labels =con$seq)
axis(side=4, at=c(1:22), labels=aa)
grid(length(con$seq), length(aa))
box()

# Add consensus sequence
for(i in 1:length(con$seq)) {
  text(i, which(aa==con$seq[i]),con$seq[i],col="white")
}

# Add lines for residue type separation
abline(h=c(2.5,3.5, 4.5, 5.5, 3.5, 7.5, 9.5,
         12.5, 14.5, 16.5, 19.5), col="gray")


}

\keyword{ utilities }
