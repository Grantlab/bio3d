\name{var.xyz}
\alias{var.xyz}
\alias{var.pdbs}
\title{ Pairwise Distance Variance in Cartesian Coordinates }
\description{
  Calculate the variance of all pairwise distances in an ensemble of
  Cartesian coordinates.
}
\usage{
var.xyz(xyz, weights=TRUE)
var.pdbs(pdbs, ...)
}
\arguments{
  \item{xyz}{ an object of class \code{"xyz"} containing Cartesian
    coordinates in a matrix. }
  \item{weights }{ logical, if TRUE weights are calculated based on the
    pairwise distance variance. }
  \item{pdbs}{ a \sQuote{pdbs} object as object from function
    \code{pdbaln}. }
  \item{\dots}{ arguments passed to associated functions. }
}
\details{
  This function calculates the variance of all pairwise distances in an
  ensemble of Cartesian coordinates. The primary use of this function is
  to calculate weights to scale the pair force constant for NMA. 
}
\value{
  Returns the a matrix of the pairwise distance variance, formated as
  weights if \sQuote{weights=TRUE}. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma.pdbs}} }
\keyword{ utilities }
