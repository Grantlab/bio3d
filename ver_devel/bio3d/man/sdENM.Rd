\name{sdENM}
\alias{sdENM}
\title{ Index for the sdENM ff }
\description{
  A dictonary of spring force constants for the sdENM force field.
}
\usage{
data(sdENM)
}
\format{
 A data frame containing the following columns:
 \enumerate{
   \item{aa1}{ one-letter code of amino acid residue i }
   \item{aa2}{ one-letter code of amino acid residue j }
   \item{min}{ minimum distance in which the corresponding spring force
     constant is value }
   \item{max}{ maximum distance in which the corresponding spring force
     constant is value  }
   \item{k}{ spring force constant valid for the corresponding distance
     range and amino acid pair. } 
}
}
\source{
  Dehouck Y. & Mikhailov A.S. (2013) \emph{PLoS Comput Biol} \bold{9}:e1003209.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\examples{
## Load force constant data
data(sdENM)
}
\keyword{datasets}
