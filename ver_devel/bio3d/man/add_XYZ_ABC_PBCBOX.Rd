\encoding{UTF-8}
\name{addXYZ}
\alias{addXYZ}
\alias{addABC}
\alias{addPBCBox}
\title{Add Cartesian Axes, Lattice Vectors or a PBC box to a Scene}
\description{Add Cartesian axes, lattice vectors or a PBC box to the current scene}
\usage{
addXYZ(lwd = 2, labels = TRUE, cex = 2)

addABC(x, lwd = 2, labels = TRUE, cex = 2)

addPBCBox(x, lwd = 2)
}
\arguments{
  \item{lwd}{a numeric value indicating the line width used to draw the Cartesian axes, lattice vectors and PBC box.}
  \item{labels}{a logical value indicating whether the labels of the Cartesian axes and lattice vectors have to drawn.}
  \item{cex}{a numeric value indicating the magnification used to draw the labels of the Cartesian axes and lattice vectors.}
  \item{x}{a length 6 numeric vector containing the a, b, c, alpha, beta and gamma lattice parameters.}
}
\details{
addXYZ:\cr
Add Cartesian axes x, y and z to the current rgl device.\cr
addABC:\cr
Add lattice vectors a, b and c to the current rgl device.\cr
addPBCBox:\cr
Add a box representing the Periodic Boundary Conditions of a molecular system.
}
\value{
Return (using invisible) a two-column data.frame containing the IDs
and type indicators of the objects added to the scene.
}
\author{Julien Idé}

\seealso{
\code{\link{visualize}}, \code{rgl.cur}, \code{par3d}
}
\examples{
## Read a PDB file from the RCSB online database
pdb <- read.pdb("1bg2")
visualize(pdb, xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE)
## Add Cartesian axes
addXYZ()
## Add cubic lattice vectors
cell <- c(70,70,70,90,90,90)
addABC(cell)
## Add a cubic PBC box
addPBCBox(cell)
}
\keyword{ dynamic }
