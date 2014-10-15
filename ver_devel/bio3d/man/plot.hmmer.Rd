\name{plot.hmmer}
\alias{plot.hmmer}
\title{ Plot a Summary of HMMER Hit Statistics. }
\description{
  Produces a number of basic plots that should facilitate hit selection
  from the match statistics of a HMMER result.
}
\usage{
 \method{plot}{hmmer}(x, cutoff = NULL, cut.seed=NULL, cluster=TRUE, mar=c(2, 5, 1, 1), cex=1.1, ...)
}
\arguments{
  \item{x}{ HMMER results as obtained from the function
    \code{\link{hmmer}}. }
  \item{cutoff}{ A numeric cutoff value, in terms of minus the log of
    the evalue, for returned hits.  If null then the function will
  try to find a suitable cutoff near \sQuote{cut.seed} which can be used
  as an initial guide (see below). }
  \item{cut.seed}{ A numeric seed cutoff value, used for initial
    cutoff estimation. If null then a seed position is set to the point of 
    largest drop-off in normalized scores (i.e. the biggest jump in 
    E-values). }
  \item{cluster}{ Logical, if TRUE (and \sQuote{cutoff} is null) a clustering 
    of normalized scores is performed to partition hits in groups by 
    similarity to query. If FALSE the partition point is set to the point of largest drop-off in normalized scores. }
  \item{mar}{ A numerical vector of the form c(bottom, left, top, right)
    which gives the number of lines of margin to be specified on the
    four sides of the plot.}
  \item{cex}{ a numerical single element vector giving the amount by
    which plot labels should be magnified relative to the default. }
  \item{\dots}{ extra plotting arguments. }  
}
\details{
  Examining plots of HMMER scores, E-values and normalized
  scores (-log(E-Value), see \sQuote{hmmer} function) can aid in the
  identification sensible hit similarity thresholds.

  If a \sQuote{cutoff} value is not supplied then a basic hierarchical
  clustering of normalized scores is performed with initial group partitioning
  implemented at a hopefully sensible point in the vicinity of
  \sQuote{h=cut.seed}. Inspection of the resultant plot can then be use
  to refine the value of \sQuote{cut.seed} or indeed
  \sQuote{cutoff}. As the \sQuote{cutoff} value can vary depending on
  the desired application and indeed the properties of the system under
  study it is envisaged that \sQuote{plot.hmmer} will be called multiple
  times to aid selection of a suitable \sQuote{cutoff} value. See the
  examples below for further details.
}
\value{
  Produces a plot on the active graphics device and returns a three
  component list object:
  \item{hits}{ an ordered matrix detailing the subset of hits with a normalized
    score above the chosen cutoff. Database identifiers are listed along
    with their cluster group number. }
  \item{acc}{ a character vector containing the database accession identifier of
    each hit above the chosen threshold. }
  \item{inds}{ a numeric vector containing the indices of the hits
               relative to the input hmmer object.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\note{ TO BE IMPROVED. 
}
\seealso{ \code{\link{hmmer}}, \code{\link{blast.pdb}} }
\examples{
\dontrun{
##- PHMMER
seq <- get.seq("2abl_A", outfile = tempfile())
res <- hmmer(seq, db="pdb")
plot.hmmer(res)
}
}
\keyword{ hplot }
