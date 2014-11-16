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
  An array of 27 matrices containg the spring force constants for the
  \sQuote{sdENM} force field (see Dehouch et al for more
  information). Each matrix in the array holds the force constants for
  all amino acid pairs for a specific distance range.

  See examples for more details.
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

## force constants for amino acids A, C, D, E, and F
## in distance range [4, 4.5)
sdENM[1:5, 1:5, 1]

## and distance range [4.5, 5)
sdENM[1:5, 1:5, 2]

## amino acid pair A-P, at distance 4.2
sdENM["A", "P", 1]

\dontrun{
## for use in NMA
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
modes <- nma(pdb, ff="sdenm")
}
}
\keyword{datasets}
