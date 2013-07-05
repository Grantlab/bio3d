\name{atom.index}
\alias{atom.index}
\title{ Index of Atomic Masses }
\description{
  A dictonary of atomic masses and PDB atom names to atom elements
  mapping.
}
\usage{
data(atom.index)
}
\format{
 A list with the following components:
 \enumerate{
  \item{elety}{ list: PDB atom name indices and associated atom
    element types. }
  \item{mass}{ list: relative atomic masses. }
}
}
\source{
  Most text books in chemistry.\cr
  See also: \url{http://en.wikipedia.org/wiki/Relative_atomic_mass}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\examples{
## Load atom index data
data(atom.index)

# Examine the masses
atom.index$mass

# Fetch mass for an atom element
atom.index$mass[["C"]]

# Examine atom name to atom element mappings
atom.index$elety

# Fetch atom element for a specific PDB atom name
atom.index$elety[["SD"]]

}
\keyword{datasets}
