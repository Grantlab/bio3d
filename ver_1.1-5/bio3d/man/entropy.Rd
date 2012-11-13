\name{entropy}
\alias{entropy}
\title{ Shannon Entropy Score }
\description{
  Calculate the sequence entropy score for every position in an alignment.
}
\usage{
entropy(alignment)
}
\arguments{
  \item{alignment}{ sequence alignment returned from
    \code{\link{read.fasta}} or an alignment character matrix. }
}
\details{
  Shannon's information theoretic entropy (Shannon, 1948) is an
  often-used measure of residue diversity and hence residue conservation.
}
\value{
  Returns a list with five components:
  \item{H }{standard entropy score for a 22-letter alphabet.}
  \item{H.10 }{entropy score for a 10-letter alphabet (see below).}
  \item{H.norm }{ normalized entropy score (for 22-letter alphabet), so that 
    conserved (low entropy) columns (or positions) score 1, and diverse
    (high entropy) columns score 0.}
  \item{H.10.norm }{ normalized entropy score (for 10-letter alphabet),
    so that conserved (low entropy) columns score 1 and diverse
    (high entropy) columns score 0.}  
  \item{freq }{residue frequency matrix containing percent occurrence
    values for each residue type.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Shannon (1948) \emph{The System Technical J.} \bold{27}, 379--422.
  
  Mirny and Shakhnovich (1999) \emph{J. Mol. Biol.} \bold{291}, 177--196.
}
\author{ Barry Grant }
\note{ 
  In addition to the standard entropy score (based on a 22-letter
  alphabet of the 20 standard amino-acids, plus a gap character \sQuote{-}
  and a mask character \sQuote{X}), an entropy score, \code{H.10}, based on
  a 10-letter alphabet is also returned.

  For \code{H.10}, residues from the 22-letter alphabet are classified
  into one of 10 types, loosely following the convention of Mirny and
  Shakhnovich (1999):
  Hydrophobic/Aliphatic [V,I,L,M],
  Aromatic [F,W,Y],
  Ser/Thr [S,T],
  Polar [N,Q],
  Positive [H,K,R],
  Negative [D,E],
  Tiny [A,G],
  Proline [P],
  Cysteine [C], and
  Gaps [-,X].

  The residue code \sQuote{X} is useful for handling non-standard aminoacids.
}
\seealso{ \code{\link{consensus}}, \code{\link{read.fasta}} }
\examples{

# Read HIV protease alignment 
aln <- read.fasta(system.file("examples/hivp_xray.fa",package="bio3d"))

# Entropy and consensus
h   <- entropy(aln)
con <- consensus(aln)

names(h$H)=con$seq
print(h$H)

# Entropy for sub-alignment (positions 1 to 20) 
h.sub <- entropy(aln$ali[,1:20])

# Plot entropy and residue frequencies (excluding positions >=60 percent gaps)
H <- h$H.norm
H[ apply(h$freq[21:22,],2,sum)>=0.6 ] = 0

col <- mono.colors(32)
aa  <- rev(rownames(h$freq))
oldpar <- par(no.readonly=TRUE)
layout(matrix(c(1,2),2,1,byrow = TRUE), widths = 7, 
       heights = c(2, 8), respect = FALSE)

# Plot 1: entropy
par(mar = c(0, 4, 2, 2))
barplot(H, border="white", ylab = "Entropy",
        space=0, xlim=c(3.7, 97.3),yaxt="n" )
axis(side=2, at=c(0.2,0.4, 0.6, 0.8))
axis(side=3, at=(seq(0,length(con$seq),by=5)-0.5),
     labels=seq(0,length(con$seq),by=5))
box()

# Plot2: residue frequencies 
par(mar = c(5, 4, 0, 2))
image(x=1:ncol(con$freq),
      y=1:nrow(con$freq),
      z=as.matrix(rev(as.data.frame(t(con$freq)))),
      col=col, yaxt="n", xaxt="n",
      xlab="Alignment Position", ylab="Residue Type")
axis(side=1, at=seq(0,length(con$seq),by=5))
axis(side=2, at=c(1:22), labels=aa)
axis(side=3, at=c(1:length(con$seq)), labels =con$seq)
axis(side=4, at=c(1:22), labels=aa)
grid(length(con$seq), length(aa))
box()

for(i in 1:length(con$seq)) {
  text(i, which(aa==con$seq[i]),con$seq[i],col="white")
}
abline(h=c(3.5, 4.5, 5.5, 3.5, 7.5, 9.5,
         12.5, 14.5, 16.5, 19.5), col="gray")

par(oldpar)
}
\keyword{ utilities }
