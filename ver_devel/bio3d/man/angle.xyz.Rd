\name{angle.xyz}
\alias{angle.xyz}
\title{ Calculate the Angle Between Three Atoms }
\description{
  A function for basic bond angle determination.
}
\usage{
angle.xyz(xyz, atm.inc = 3)
}
\arguments{
  \item{xyz}{ a numeric vector of Cartisean coordinates. }
  \item{atm.inc}{ a numeric value indicating the number of atoms to
    increment by between successive angle evaluations (see below). }
}
\value{
  Returns a numeric vector of angles.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  With \code{atm.inc=1}, angles are calculated for each set of
  three successive atoms contained in \code{xyz} (i.e. moving along one
  atom, or three elements of \code{xyz}, between sucessive
  evaluations).  With \code{atm.inc=3}, angles are calculated for each set
  of three successive non-overlapping atoms contained in \code{xyz}
  (i.e. moving along three atoms, or nine elements of \code{xyz}, between
  sucessive evaluations).
}
\seealso{ \code{\link{torsion.pdb}}, \code{\link{torsion.xyz}},
  \code{\link{read.pdb}}, \code{\link{read.dcd}}.
}
\examples{
\donttest{
## Read a PDB file
pdb <- read.pdb( "http://www.rcsb.org/pdb/files/1BG2.pdb" )

## Angle between N-CA-C atoms of residue four
inds <- atom.select(pdb,"///4///N,CA,C/")
angle.xyz(pdb$xyz[inds$xyz])

## Basic stats of all N-CA-C bound angles
inds <- atom.select(pdb,"//////N,CA,C/")
summary( angle.xyz(pdb$xyz[inds$xyz]) )
#hist( angle.xyz(pdb$xyz[inds$xyz]), xlab="Angle" )

}
}
\keyword{ utilities }
