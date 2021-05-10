\name{aa2index}
\alias{aa2index}
\title{ Convert an Aminoacid Sequence to AAIndex Values }
\description{
  Converts sequences to aminoacid indeces from the \sQuote{AAindex} database. 
}
\usage{
aa2index(aa, index = "KYTJ820101", window = 1)
}
\arguments{
  \item{aa}{ a protein sequence character vector. }
  \item{index}{ an index name or number (default: \dQuote{KYTJ820101},
    hydropathy index by Kyte-Doolittle, 1982). }
  \item{window}{ a positive numeric value, indicating the window size
    for smoothing with a sliding window average (default: 1, i.e. no
    smoothing). }
}
\details{
  By default, this function simply returns the index values for each
  amino acid in the sequence. It can also be set to perform a crude
  sliding window average through the \code{window} argument. 
}
\value{
  Returns a numeric vector.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  \sQuote{AAIndex} is the work of Kanehisa and co-workers:
Kawashima and Kanehisa (2000) \emph{Nucleic Acids Res.} \bold{28}, 374;
Tomii and Kanehisa (1996) \emph{Protein Eng.} \bold{9}, 27--36;
Nakai, Kidera and Kanehisa (1988) \emph{Protein Eng.} \bold{2}, 93--100.

For a description of the \sQuote{AAindex} database see:\cr 
\url{https://www.genome.jp/aaindex/} or the \code{\link{aa.index}} documentation.
}
\author{ Ana Rodrigues }
\seealso{ \code{\link{aa.index}}, \code{\link{read.fasta}} }
\examples{
## Residue hydropathy values
seq <- c("R","S","D","X","-","X","R","H","Q","V","L")
aa2index(seq)

\dontrun{
## Use a sliding window average
aa2index(aa=seq, index=22, window=3)

## Use an alignment

aln  <- read.fasta(system.file("examples/hivp_xray.fa",package="bio3d"))
prop <- t(apply(aln$ali, 1, aa2index, window=1))

## find and use indices for volume calculations
i <- which(sapply(aa.index,
       function(x) length(grep("volume", x$D, ignore.case=TRUE)) != 0))
sapply(i, function(x) aa2index(aa=seq, index=x, window=5)) 
}
}
\keyword{ utilities }

