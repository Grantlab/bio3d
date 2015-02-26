\name{aa.table}
\alias{aa.table}
\docType{data}
\title{ Table of Relevant Amino Acids }
\description{
  This data set provides the atomic masses of a selection of amino
  acids regularly occuring in proteins. 
}
\usage{
aa.table
}
\format{
  A data frame with the following components.
  \describe{
    \item{\code{aa3}}{a character vector containing three-letter amino
      acid code.}
    \item{\code{aa1}}{a character vector containing one-letter amino
      acid code.}
    \item{\code{mass}}{a numeric vector containing the mass of the
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
data(aa.table)
aa.table

## table look up
aa.table["HIS", ]

## read PDB, and fetch residue masses
pdb <- read.pdb(system.file("examples/1hel.pdb", package="bio3d"))
aa2mass(pdb)
}
\keyword{datasets}
