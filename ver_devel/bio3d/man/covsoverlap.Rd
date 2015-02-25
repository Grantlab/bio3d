\name{covsoverlap}
\alias{covsoverlap}
\alias{covsoverlap.enma}
\alias{covsoverlap.nma}
\title{ Covariance Overlap }
\description{
  Calculate the covariance overlap obtained from NMA. 
}
\usage{
covsoverlap(...)

\method{covsoverlap}{enma}(enma, ncore=NULL, ...)

\method{covsoverlap}{nma}(a, b, subset=NULL, ...)
}
\arguments{
  \item{enma}{ an object of class \code{"enma"} obtained from function
    \code{nma.pdbs}. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. } 
  \item{a}{ a list object with elements \sQuote{U} and \sQuote{L}
    (e.g. as obtained from function \code{nma}) containing the
    eigenvectors and eigenvalues, respectively, to be 
    compared with \code{b}. }
  \item{b}{ a list object with elements \sQuote{U} and \sQuote{L}
    (e.g. as obtained from function \code{nma}) containing the
    eigenvectors and eigenvalues, respectively, to be 
    compared with \code{a}. }
  \item{subset}{ the number of modes to consider. }
  \item{\dots}{ arguments passed to associated functions. }
}
\details{
  Covariance overlap is a measure for the similarity between two
  covariance matrices, e.g. obtained from NMA. 
}
\value{
  Returns the similarity coefficient(s).
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
  Romo, T.D. et al. (2011) \emph{Proteins} \bold{79}, 23--34.
}
\author{ Lars Skjaerven }
\seealso{ Other similarity measures:
  \code{\link{sip}}, \code{\link{covsoverlap}}, \code{\link{bhattacharyya}}.
}
\keyword{ utilities }
