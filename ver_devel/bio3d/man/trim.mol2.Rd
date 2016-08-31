\name{trim.mol2}
\alias{trim.mol2}
\title{ Trim a MOL2 Object To A Subset of Atoms. }
\description{
  Produce a new smaller MOL2 object, containing a subset of atoms, from a
  given larger MOL2 object. 
}
\usage{
\method{trim}{mol2}(mol, \dots, inds = NULL)
}
\arguments{
  \item{mol}{ a MOL2 structure object obtained from
    \code{\link{read.mol2}}. }
  \item{\dots}{ additional arguments passed to \code{\link{atom.select}}. If
     \code{inds} is also provided, these arguments will be ignored. }
  \item{inds}{ a list object of ATOM and XYZ indices as obtained from
    \code{\link{atom.select}}. If NULL, atom selection will be obtained
    from calling \code{atom.select(mol, \dots)}. }
}
\details{
  This is a basic utility function for creating a new MOL2 object based
  on a selection of atoms.
}
\value{
  Returns a list of class \code{"mol2"}.
  }
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{read.mol2}}, \code{\link{atom.select.mol2}},
  \code{\link{as.pdb.mol2}}
}
\examples{
\dontrun{
## Read a MOL2 file from those included with the package
mol <- read.mol2( system.file("examples/aspirin.mol2", package="bio3d"))

## Trim away H-atoms
mol <- trim(mol, "noh")

}
}
\keyword{ utilities }
