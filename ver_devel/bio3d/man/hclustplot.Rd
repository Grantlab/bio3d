\name{hclustplot}
\alias{hclustplot}
\title{ Dendrogram with Clustering Annotation }
\description{
  Draw a standard dendrogram with clustering annotation in the marginal
  regions and colored labels.
}
\usage{
hclustplot(hc, k = NULL, h = NULL, colors = NULL, labels = NULL,
             fillbox = TRUE, heights = c(1, .3), mar = c(1, 1, 0, 1), ...)
}
\arguments{
  \item{hc}{ an object of the type produced by \code{hclust}. }
  \item{k}{ an integer scalar or vector with the desired number of
    groups. Redirected to function \code{cutree}. }
  \item{h}{ numeric scalar or vector with heights where the tree should
    be cut. Redirected to function \code{cutree}. At least one of \sQuote{k}
    or \sQuote{h} must be specified. }
  \item{colors}{ a numerical or character vector with the same length as
    \sQuote{hc} specifying the colors of the labels. }
  \item{labels}{ a character vector with the same length as
    \sQuote{hc} containing the labels to be written. }
  \item{fillbox}{ logical, if TRUE clustering annotation will be drawn
    as filled boxes below the dendrogram. }
  \item{heights}{ numeric vector of length two specifying the values for
    the heights of rows on the device. See function \code{layout}. }
  \item{mar}{ a numerical vector of the form \sQuote{c(bottom, left, top, right)}
    which gives the number of lines of margin to be specified on
    the four sides of the plot. If left at default the margins will be
    adjusted upon adding arguments \sQuote{main}, \sQuote{ylab}, etc. }
  \item{\dots}{ other graphical parameters passed to functions
    \code{plot.dendrogram}, \code{mtext}, and \code{par}. Note that
    certain arguments will be ignored. }
}
\details{
  This function adds extended visualization of cluster membership to a
  standard dendrogram. If \sQuote{k} or \sQuote{h} is provided a call to
  \code{cutree} will provide cluster membership
  information. Alternatively a vector of colors or cluster membership
  information can be provided through argument \sQuote{colors}.

  See examples for further details on usage.
}
\note{
  Argument \sQuote{horiz=TRUE} currently not supported.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{plot.hclust}}, \code{\link{plot.dendrogram}},
  \code{\link{hclust}}, \code{\link{cutree}}.  }
\examples{

attach(transducin)

##- perform clustering
rd <- rmsd(pdbs, fit=TRUE)
hc <- hclust(as.dist(rd))

##- draw dendrogram 
hclustplot(hc, k=3)

##- draw dendrogram with manuall clustering annotation
hclustplot(hc, colors=annotation[, "color"], labels=pdbs$id)

##- draw dendrogram 
hclustplot(hc, k=3, labels=pdbs$id, fillbox=FALSE)



}
\keyword{ hplot }
