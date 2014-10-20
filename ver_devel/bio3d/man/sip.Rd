\name{sip}
\alias{sip}
\alias{sip.default}
\alias{sip.nma}
\alias{sip.enma}
\title{ Square Inner Product }
\description{
  Calculate the correlation between two atomic fluctuation vectors.
}
\usage{
sip(...)
\method{sip}{nma}(a, b, ...)
\method{sip}{enma}(enma, ncore=NULL, ...)
\method{sip}{default}(v, w, ...)
}
\arguments{
  \item{enma}{ an object of class \code{"enma"} obtained from function
    \code{nma.pdbs}. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{a}{ an \sQuote{nma} object as object from function \code{nma} to
    be compared to \code{b}. }
  \item{b}{ an \sQuote{nma} object as object from function \code{nma} to
    be compared to \code{a}. }
  \item{v}{ a numeric vector containing the atomic fluctuation values. }
  \item{w}{ a numeric vector containing the atomic fluctuation values. }
  \item{\dots}{ arguments passed to associated functions. }
}
\details{
  SIP is a measure for the similarity of atomic fluctuations of two
  proteins, e.g. experimental b-factors, theroetical RMSF values, or 
  atomic fluctuations obtained from NMA.
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
  \code{\link{covsoverlap}}, \code{\link{bhattacharyya}},
  \code{\link{rmsip}}.
}
\examples{

pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
a <- nma(pdb)
b <- nma(pdb, ff="anm")

sip(a$fluctuations, b$fluctuations)

}
\keyword{ utilities }
