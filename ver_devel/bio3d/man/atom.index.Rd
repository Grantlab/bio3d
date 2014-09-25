\name{atom.index}
\alias{atom.index}
\docType{data}
\title{
  Atom Names/Types
}
\description{
  This data set gives for various atom names/types the corresponding
  atomic symbols.
}
\usage{
atom.index
}
\format{
  A data frame with the following components.
  \describe{
    \item{\code{name}}{a character vector containing atom names/types.}
    \item{\code{symb}}{a character vector containing atomic symbols.}
  }
}
\seealso{
  \code{\link{elements}}, \code{\link{atom.index}},
  \code{\link{atom2ele}}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\examples{
data(atom.index)
atom.index

# Get the atomic symbol of some atoms
atom.names <- c("CA", "O", "N", "OXT")
atom.index[match(atom.names, atom.index$name), "symb"]

}
\keyword{datasets}
