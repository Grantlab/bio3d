\name{plot.hmmer}
\alias{plot.hmmer}
\title{ Plot a Summary of HMMER Hit Statistics. }
\description{
  Produces a number of basic plots that should facilitate hit selection
  from the match statistics of a HMMER result.
}
\usage{
 \method{plot}{hmmer}(x, ...)
}
\arguments{
  \item{x}{ HMMER results as obtained from the function
    \code{\link{hmmer}}. }
  \item{\dots}{ arguments passed to \code{\link{plot.blast}}. }
}
\details{
  See \code{\link{plot.blast}} for details. 
}
\value{
  Produces a plot on the active graphics device and returns a three
  component list object:
  \item{hits}{ an ordered matrix detailing the subset of hits with a normalized
    score above the chosen cutoff. Database identifiers are listed along
    with their cluster group number. }
  \item{acc}{ a character vector containing the database accession identifier of
    each hit above the chosen threshold. }
  \item{pdb.id}{ a character vector containing the database accession identifier of
    each hit above the chosen threshold. }
  \item{inds}{ a numeric vector containing the indices of the hits
               relative to the input hmmer object.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{ \code{\link{hmmer}}, \code{\link{blast.pdb}} }
\examples{
\dontrun{
# HMMER server connection required - testing excluded

##- PHMMER
seq <- get.seq("2abl_A", outfile = tempfile())
res <- hmmer(seq, db="pdb")
plot.hmmer(res)

}
}
\keyword{ hplot }
