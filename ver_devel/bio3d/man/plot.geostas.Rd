\name{plot.geostas}
\alias{plot.geostas}
\title{ Plot Geostas Results }
\description{
  Plot an atomic movement similarity matrix with domain annotation
}
\usage{
\method{plot}{geostas}(x, at=seq(0, 1, 0.1), main="AMSM with Domain Assignment",
                         col.regions=rev(heat.colors(200)), 
                         margin.segments=x$grps, ...)
}
\arguments{
  \item{x}{ an object of type \code{geostas} as obtained by the
    \sQuote{geostas} function. } 
  \item{at}{ numeric vector specifying the levels to be colored. }
  \item{main}{ a main title for the plot. }
  \item{col.regions}{ color vector. See \code{contourplot} for more 
    information. }
  \item{margin.segments}{ a numeric vector of cluster membership as 
    obtained from cutree() or other community detection method. This 
    will be used for bottom and left margin annotation. }
  \item{\dots}{ additional graphical parameters for
    \code{\link{plot.dccm}} and \code{contourplot}. }
}
\details{
  This is a wrapper function for \code{\link{plot.dccm}} with appropriate
  adjustments for plotting atomic movement similarity matrix obtained
  from function \code{\link{geostas}}. 
  
  See the \code{\link{plot.dccm}} for more details.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{
  \code{\link{plot.dccm}}, \code{\link{geostas}}
}
\keyword{ hplot }
