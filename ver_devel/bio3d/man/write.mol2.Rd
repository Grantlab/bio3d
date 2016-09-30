\name{write.mol2}
\alias{write.mol2}
\title{ Write MOL2 Format Coordinate File }
\description{
  Write a Sybyl MOL2 file
}
\usage{
write.mol2(mol, file = "R.mol2",  append = FALSE)
}
\arguments{
  \item{mol}{ a MOL2 structure object obtained from
    \code{\link{read.mol2}}. }
  \item{file}{ the output file name. }
  \item{append}{ logical, if TRUE output is appended to the bottom of an
    existing file (used primarly for writing multi-model files). }
}
\details{
  See examples for further details.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{
  Lars Skjaerven
}
\examples{
\donttest{
# Read MOL2 file
mol <- read.mol2( system.file("examples/aspirin.mol2", package="bio3d") )

# Trim away H-atoms
mol <- trim(mol, "noh")

# Write new MOL2 file
#write.mol2(mol)
}
}
\keyword{ IO }
