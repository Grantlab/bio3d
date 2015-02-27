\name{aa.index}
\alias{aa.index}
\docType{data}
\title{ AAindex: Amino Acid Index Database }
\description{
  A collection of published indices, or scales, of numerous
  physicochemical and biological properties of the 20 standard
  aminoacids (Release 9.1, August 2006).
}
\usage{data(aa.index)}
\format{
 A list of 544 named indeces each with the following components:
 \enumerate{
  \item{H}{ character vector: Accession number. }
  \item{D}{ character vector: Data description. }
  \item{R}{ character vector: LITDB entry number. }
  \item{A}{ character vector: Author(s). }
  \item{T}{ character vector: Title of the article. }
  \item{J}{ character vector: Journal reference. }
  \item{C}{ named numeric vector: Correlation coefficients of similar
    indeces (with coefficients of 0.8/-0.8 or more/less). The
    correlation coefficient is calculated with zeros filled for missing
    values. }
  \item{I}{ named numeric vector: Amino acid index data. }
 }
}
\source{
\sQuote{AAIndex} was obtained from:\cr
\url{http://www.genome.jp/aaindex/}\cr
For a description of the \sQuote{AAindex} database see:\cr 
\url{http://www.genome.jp/aaindex/aaindex_help.html}.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  
\sQuote{AAIndex} is the work of Kanehisa and co-workers:\cr
Kawashima and Kanehisa (2000) \emph{Nucleic Acids Res.} \bold{28}, 374;\cr
Tomii and Kanehisa (1996) \emph{Protein Eng.} \bold{9}, 27--36;\cr
Nakai, Kidera and Kanehisa (1988) \emph{Protein Eng.} \bold{2}, 93--100.
}
\examples{
## Load AAindex data
data(aa.index)

## Find all indeces described as "volume"
ind <- which(sapply(aa.index, function(x)
                    length(grep("volume", x$D, ignore.case=TRUE)) != 0))

## find all indeces with author "Kyte"
ind <- which(sapply(aa.index, function(x) length(grep("Kyte", x$A)) != 0))

## examine the index
aa.index[[ind]]$I

## find indeces which correlate with it
all.ind <- names(which(Mod(aa.index[[ind]]$C) >= 0.88))

## examine them all
sapply(all.ind, function (x) aa.index[[x]]$I)
}
\keyword{datasets}
