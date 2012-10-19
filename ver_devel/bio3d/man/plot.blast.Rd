\name{plot.blast}
\alias{plot.blast}
\title{ Plot a Summary of BLAST Hit Statistics. }
\description{
  Produces a number of basic plots that should facilitate hit selection
  from the match statistics of a BLAST result.
}
\usage{
plot.blast(x, cutoff = NULL, cut.seed=110, mar=c(4, 4, 1, 2), cex.lab=1.5, ...)
}
\arguments{
  \item{x}{ BLAST results as obtained from the function
    \code{\link{blast.pdb}}. }
  \item{cutoff}{ A numeric cutoff value, in terms of minus the log of
    the evalue, for returned hits.  If null then the function will
  try to find a suitable cutoff near \sQuote{cut.seed} which can be used
  as an initial guide (see below). }
  \item{cut.seed}{ A numeric seed cutoff value, used for initial
    cutoff estimation. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{cex.lab}{ a numerical single element vector giving the amount by
    which plot labels should be magnified relative to the default. }
  \item{\dots}{ extra plotting arguments. }  
}
\details{
  Examining plots of BLAST alignment lengths, scores, E-values and normalized
  scores (-log(E-Value), see \sQuote{blast.pdb} function) can aid in the
  identification sensible hit similarity thresholds.

  If a \sQuote{cutoff} value is not supplied then a basic hierarchical
  clustering of normalized scores is performed with initial group partitioning
  implemented at a hopefully sensible point in the vicinity of
  \sQuote{h=cut.seed}. Inspection of the resultant plot can then be use
  to refine the value of \sQuote{cut.seed} or indeed
  \sQuote{cutoff}. As the \sQuote{cutoff} value can vary depending on
  the desired application and indeed the properties of the system under
  study it is envisaged that \sQuote{plot.blast} will be called multiple
  times to aid selection of a suitable \sQuote{cutoff} value. See the
  examples below for further details.
}
\value{
  Produces a plot on the active graphics device and returns a three
  component list object:
  \item{hits}{ an ordered matrix detailing the subset of hits with a normalized
    score above the chosen cutoff. Database identifiers are listed along
    with their cluster group number. }
  \item{pdb.id}{ a character vector containing the PDB database identifier of
    each hit above the chosen threshold. }
  \item{gi.id}{ a character vector containing the gi database identifier of
    each hit above the chosen threshold. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ TO BE IMPROVED. 
}
\seealso{ \code{\link{blast.pdb}} }
\examples{
b2 <-  blast.pdb( seq.pdb(read.pdb( get.pdb("4q21", URLonly=TRUE) )) )
raw.hits <- plot.blast(b2)
top.hits <- plot.blast(b2, 188)
head(top.hits$hits)

\dontrun{
blast <- blast.pdb( seq.pdb(read.pdb( get.pdb("2BN3", URLonly=TRUE) )))
raw.hits <- plot(blast)
top.hits <- plot(blast, cut.seed=20)

head(top.hits$pdb.id)
pdbFiles <- get.pdb(top.hits$pdb.id, dir="downloadedPDBs/")
split.pdb(pdbFiles, dir="PDB_chains/")
}
}
\keyword{ hplot }
