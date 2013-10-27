\name{dccm}
\alias{dccm}
\title{ DCCM: Dynamical Cross-Correlation Matrix }
\description{
  Determine the cross-correlations of atomic displacements. 
}
\usage{
dccm(x, ...)
}
\arguments{
  \item{x}{ a numeric matrix of Cartesian coordinates with a row per
    structure/frame which will br passed to \code{dccm.xyz()}.
    Alternatively, an object  of class \code{nma} as obtained from
    function \code{nma} that will be passed to the \code{dccm.nma()}
    function, see below for examples. } 
  \item{\dots}{ additional arguments passed to the methods
    \code{dccm.xyz}, \code{dccm.nma}, and \code{dccm.enma}. }
}
\details{
  \sQuote{dccm} is a generic function...
}
\value{
  Returns a cross-correlation matrix.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{ \code{\link{dccm.xyz}}, \code{\link{dccm.nma}},
  \code{\link{dccm.enma}} . }
\keyword{ utilities }
