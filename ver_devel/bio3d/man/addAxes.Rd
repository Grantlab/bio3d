\name{addAxes}
\alias{addABC}
\alias{addXYZ}
\alias{addPBCBox}
\title{Add Axes or PBC Box to the \sQuote{rgl} Scene}
\description{Add lattice vectors, Cartesian axes or PBC box to the current \sQuote{rgl} scene.}
\usage{
  addABC(x, lwd = 2, labels = TRUE, cex = 2)
  addXYZ(lwd = 2, labels = TRUE, cex = 2, col = "white")
  addPBCBox(x, lwd = 2, col = "white")
}
\arguments{
  \item{x}{an object of class \sQuote{cryst1} containing unit cell parameters.}
  \item{lwd}{a numeric value indicating the line width used to draw the axes or the PBC box.}
  \item{labels}{a logical value indicating whether the labels of the axes have to be drawn.}
  \item{cex}{a numeric value indicating the magnification used to draw the labels of the axes.}
  \item{col}{the color used to draw Cartesian the axes or the PBC box.}
}
\details{
  \code{addABC}: Add the lattice vectors a, b and c to the current rgl device.\cr
  \code{addXYZ}: Add the Cartesian axes x, y and z to the current rgl device.\cr
  \code{addPBCBox}: Add a box representing the Periodic Boundary Conditions of a molecular system.
}
\value{
  Return (using invisible) a two-column data.frame containing the IDs
  and type indicators of the objects added to the scene.
}
\seealso{
  \code{\link{visualize}}
}
\examples{
\dontrun{
  x <- read.pdb("1bg2")
  visualize(x, type = "l", xyz.axes = , abc.axes = FALSE, pbc.box = FALSE)
  addXYZ()
  addABC(x$cell)
  addPBCBox(x$cell)
}
}
\keyword{ dynamic }
