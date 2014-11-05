\name{aa.mass}
\alias{aa.mass}
\docType{data}
\title{ Amino Acid Residue Mass }
\description{
  This data set provides the atomic masses of a selection of amino
  acids regularly occuring in proteins. 
}
\usage{
aa.mass
}
\format{
  A data frame with the following components.
  \describe{
    \item{\code{aa3}}{a character vector containing three-letter amino
      acid code.}
    \item{\code{aa1}}{a character vector containing one-letter amino
      acid code.}
    \item{\code{aaMass}}{a numeric vector containing the mass of the
      respective amino acids. }
    \item{\code{formula}}{a character vector containing the formula of
      the amino acid in which the mass calculat was based. }
    \item{\code{name}}{a character vector containing the full names of
      the respective amino acids. }
  }
}
\seealso{
  \code{\link{aa2mass}}, \code{\link{aa.index}},
  \code{\link{atom.index}}, \code{\link{elements}}, 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\examples{
\donttest{
data(aa.mass)
aa.mass

## table look up
aa.mass["HIS", ]

## read PDB, and fetch residue masses
pdb <- read.pdb("1etl")
aa2mass(pdb)
}
}
\keyword{datasets}
