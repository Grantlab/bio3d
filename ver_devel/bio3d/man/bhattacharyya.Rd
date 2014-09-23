\name{bhattacharyya}
\alias{bhattacharyya}
\alias{bhattacharyya.nma}
\alias{bhattacharyya.pca}
\alias{bhattacharyya.enma}
\alias{bhattacharyya.array}
\alias{bhattacharyya.matrix}
\title{ Bhattacharyya Coefficient }
\description{
  Calculate the Bhattacharyya Coefficient as a similarity between two
  modes objects.
}
\usage{
bhattacharyya(...)
\method{bhattacharyya}{enma}(enma, covs=NULL, ncore=NULL, ...)
\method{bhattacharyya}{array}(covs, ncore=NULL, ...)
\method{bhattacharyya}{matrix}(a, b, q=90, n=NULL, ...)
\method{bhattacharyya}{nma}(...)
\method{bhattacharyya}{pca}(...)
}
\arguments{
  \item{enma}{ an object of class \code{"enma"} obtained from function
    \code{nma.pdbs}. }
  \item{covs}{ an array of covariance matrices of equal dimensions. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. } 
  \item{a}{ covariance matrix to be compared with \code{b}. }
  \item{b}{ covariance matrix to be compared with \code{a}. }
  \item{q}{ a numeric value (in percent) determining the number of modes
    to be compared. }
  \item{n}{ the number of modes to be compared. }
  \item{\dots}{ arguments passed to associated functions. }
}
\details{
  Bhattacharyya coefficient provides a means to compare two covariance
  matrices derived from NMA or an ensemble of conformers
  (e.g. simulation or X-ray conformers).  
}
\value{
  Returns the similarity coefficient(s).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
  Fuglebakk, E. et al. (2013) \emph{JCTC} \bold{9}, 5618--5628.
}
\author{ Lars Skjaerven }
\seealso{Other similarity measures:
  \code{\link{sip}}, \code{\link{covsoverlap}}, \code{\link{rmsip}}.
}
\keyword{ utilities }
