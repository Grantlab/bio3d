\name{wrap.tor}
\alias{wrap.tor}
\title{ Wrap Torsion Angle Data }
\description{
  Adjust angular data so that the absolute difference of any of the
  observations from its mean is not greater than 180 degrees.
}
\usage{
wrap.tor(data, wrapav=TRUE, avestruc=NULL)
}
\arguments{
  \item{data}{ a numeric vector or matrix of torsion angle data as
    obtained from \code{torsion.xyz}. }
  \item{wrapav}{ logical, if TRUE average structure is also
    \sQuote{wrapped} }
  \item{avestruc}{ a numeric vector corresponding to the average
    structure } 
  
}
\details{
  This is a basic utility function for coping with the periodicity of 
  torsion angle data, by \sQuote{wraping} angular data such that the 
  absolute difference of any of the observations from its column-wise 
  mean is not greater than 180 degrees.
}
\value{
  A numeric vector or matrix of wrapped torsion angle data.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Karim ElSawy }
\seealso{ \code{\link{torsion.xyz}} }
\keyword{ utilities }
