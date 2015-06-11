\name{cell}
\alias{cell}
\alias{cell.default}
\alias{cell.pdb}
\alias{cell<-}
\alias{cell<-.pdb}
\alias{is.cell}
\alias{cellCoords}
\alias{cellCoords.cell}
\alias{cellCoords.pdb}
\alias{cellVolume}
\alias{cellVolume.cell}
\alias{cellVolume.pdb}
\alias{cellDensity}
\alias{cellDensity.cell}
\alias{cellDensity.pdb}
\title{
  Create and Manipulate \sQuote{cell} Object
}
\description{
  Create an object of class \sQuote{cell} containing the unit lattice parameters and the name of the space group to associate with an object of class \sQuote{pdb}. Compute cell properties such as the volume or the density of a unit cell.
}
\usage{
  cell(...)
  \method{cell}{default}(abc = rep(NA_real_, 3), abg = rep(90.0, 3),
     orient = "cbzy", sgroup = "P 1", \dots)
  \method{cell}{pdb}(x, ...)
  cell(x) <- value
  \method{cell}{pdb}(x) <- value
  is.cell(x)
  cellCoords(x)
  \method{cellCoords}{cell}(x)
  \method{cellCoords}{pdb}(x)
  cellVolume(x)
  \method{cellVolume}{cell}(x)
  \method{cellVolume}{pdb}(x)
  cellDensity(x)
  \method{cellDensity}{pdb}(x)
}
\arguments{
  \item{\dots}{further arguments passed to or from other methods.}
  \item{abc}{a numeric vector of length 3 containing the norms of the lattice vectors a, b and c.}
  \item{abg}{a numeric vector of length 3 containing the angles between the lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.}
  \item{orient}{a character string indicating the orientation of the lattice vectors with respecto the Cartesian reference.}
  \item{sgroup}{a character string giving the Hermann-Mauguin symbol of the space group.}
  \item{x}{an R object containing a \sQuote{cell} component to set or get, or to be tested.}
  \item{value}{an object of class \sQuote{cell}.}
}
\details{
  \code{cell} is a generic function to create objects of class \sQuote{cell}.
  The purpose of this class is to store CRYST1 records from PDB files which contain the unit cell parameters
  and the name of the space group of a molecular system stored in a PDB file.
  The default method of the \code{cell} function create an object of class \sQuote{cell} from its different components, i.e.: \code{abc}, \code{abg}, \code{orient} and \code{sgroup}.
  
  The mothods for object of class \sQuote{pdb} can be used to get or set its \code{cell} component.
  
  \code{is.cell} tests if an object of class \sQuote{cell}.
  
  \code{cellCoords} compute the coordinates of the lattice vectors in the Cartesian reference.
  
  \code{cellVolume} compute the volume of the unit cell.
  
  \code{cellDensity} compute the density of the unit cell (in g.cm^-3).
}
\value{
  \code{cell} returns a list of class \sQuote{cell} with the following components:
  \item{abc}{a numeric vector of length 3 containing the norms of the lattice vectors a, b and c.}
  \item{abg}{a numeric vector of length 3 containing the angles between the lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.}
  \item{orient}{a character string indicating the orientation of the lattice vectors with respecto the Cartesian reference.}
  \item{sgroup}{a character string giving the Hermann-Mauguin symbol of the space group.}
  
  \code{is.cell} returns TRUE if \code{x} is an object of class \sQuote{cell} and FALSE otherwise.
  
  \code{cellCoords} returns a 3x3 matrix.
  
  \code{cellVolume} returns a numeric value.
  
  \code{cellDensity} returns a numeric value.
}

\seealso{
  \code{\link{read.pdb}}, \code{\link{visualize}}
}
\examples{
  \dontrun{
  x <- read.pdb("1bg2")
  is.cell(x$cell)
  cellCoords(x)
  cellVolume(x)
  cellDensity(x)
  }
}
\keyword{ classes }
