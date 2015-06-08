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
    \code{dccm.xyz}, \code{dccm.pca}, \code{dccm.nma}, and \code{dccm.enma}. }
}
\details{
  \code{dccm} is a generic function calling the corresponding function
  determined by the class of the input argument \code{x}. Use
  \code{methods("dccm")} to get all the methods for \code{dccm}
  generic:

  \code{\link{dccm.xyz}} will be used when \code{x} is a numeric matrix
  containing Cartesian coordinates (e.g. trajectory data).
  
  \code{\link{dccm.pca}} will calculate the cross-correlations based on
  an \code{pca} object.

  \code{\link{dccm.nma}} will calculate the cross-correlations based on
  an \code{nma} object. Similarly, \code{\link{dccm.enma}} will
  calculate the correlation matrices based on an ensemble of \code{nma}
  objects (as obtained from function \code{nma.pdbs}).

  \code{\link{plot.dccm}} and \code{\link{view.dccm}} provides
  convenient functionality to plot a correlation map, and visualize the
  correlations in the structure, respectively.

  See examples for each corresponding function for more details.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{ \code{\link{dccm.xyz}}, \code{\link{dccm.nma}},
  \code{\link{dccm.enma}}, \code{\link{dccm.pca}}, \code{\link{plot.dccm}},
    \code{\link{pymol.dccm}}. }
\keyword{ utilities }
