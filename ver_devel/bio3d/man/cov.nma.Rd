\name{cov.nma}
\alias{cov.nma}
\alias{cov.enma}
\title{ Calculate Covariance Matrix from Normal Modes }
\description{
  Calculate the covariance matrix from a normal mode object.
}
\usage{
\method{cov}{nma}(nma)
\method{cov}{enma}(enma, ncore=NULL)
}
\arguments{
  \item{nma}{ an \code{nma} object as obtained from function
    \code{nma.pdb}. }
  \item{enma}{ an \code{enma} object as obtained from function
    \code{nma.pdbs}. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. } 
}
\details{
  This function calculates the covariance matrix from a \code{nma}
  object as obtained from function \code{nma.pdb} or covariance matrices
  from a \code{enma} object as obtain from function \code{nma.pdbs}. 
}
\value{
  Returns the calculated covariance matrix (function \code{cov.nma}), or
  covariance matrices (function \code{cov.enma}). 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  Fuglebakk, E. et al. (2013) \emph{JCTC} \bold{9}, 5618--5628.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}} }
\keyword{ utilities }
