\name{fluct.nma}
\alias{fluct.nma}
\title{ NMA Fluctuations }
\description{
  Calculates the atomic fluctuations from normal modes analysis.
}
\usage{
fluct.nma(nma, mode.inds=NULL)
}
\arguments{
  \item{nma}{ a list object of class \code{"nma"} (obtained with
    \code{\link{nma}}).}
  \item{mode.inds}{ a numeric vector containing the the mode numbers in
    which the calculation should be based. }
}
\details{
  Atomic fluctuations are calculated based on the \code{nma} object. By
  default all modes are included in the calculation.

  See examples for more details.
}
\value{
  Returns a numeric vector of atomic fluctuations.
}
\references{
  Hinsen, K. et al. (2000) \emph{Chemical Physics} \bold{261}, 25--37.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma} } }
\examples{
## Fetch stucture
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## Calculate (vibrational) normal modes
modes <- nma(pdb)

## Fluctuations
f <- fluct.nma(modes)

## Fluctuations of first non-trivial mode
f <- fluct.nma(modes, mode.inds=c(7,8))
}
\keyword{ analysis }
