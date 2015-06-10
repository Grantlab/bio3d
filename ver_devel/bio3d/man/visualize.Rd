\encoding{UTF-8}
\name{visualize}
\alias{visualize}
\alias{visualize.default}
\alias{visualize.xyz}
\alias{visualize.pdb}
\alias{visualize.character}
\title{Visualize a Molecular Structure}
\description{Use the rgl library to visualize in 3D a molecular structure.}
\usage{
visualize(...)

\method{visualize}{xyz}(
xyz, elesy = NULL, con = NULL, cell = NULL, type = "l",
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rvdw", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...)

\method{visualize}{pdb}(
  pdb, cell = NULL, type = "l", con = TRUE,
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rvdw", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...)
  
\method{visualize}{character}(
  file, cell = NULL, type = "l", con = TRUE,
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rcov", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...)

}
\arguments{
  \item{pdb}{an object of class \sQuote{pdb} containing the molecular structure to visualize.}
  \item{xyz}{a numeric vector of Cartesian coordinates like the \code{xyz} component of objects of class \sQuote{pdb}.}
  \item{elesy}{a character vector containing atomic symbols (If \code{length(ele.symb)!=length(xyz)/3} \code{elesy} is recycled).}
  \item{con}{an object of class \sQuote{connectivity} or a logical value indicating if the connectivity has to be calculated}
  \item{cell}{a length 6 numeric vector containing the a, b, c, alpha, beta and gamma lattice parameters.}
  \item{type}{a character string indicating the visualization mode (See details).}
  \item{xyz.axes}{a logical value indicating whether the x, y and z axes have to be added to the scene.}
  \item{abc.axes}{a logical value indicating whether the lattice vectors have to be added to the scene.}
  \item{pbc.box}{a logical value indicating whether the pbc box has to be added to the scene.}
  \item{centre}{a logical value indicating if the structure as to be translated to its barrycentre.}
  \item{lwd}{a numeric value indication the line width used to plot the axes, the pbc box and atomic bonds when \code{type = "l"}.}
  \item{lwd.xyz}{a numeric value indicating the line width used to plot the x, y and z axes.}
  \item{lwd.abc}{a numeric value indicating the line width used to plot the lattice vectors.}
  \item{lwd.pbc.box}{a numeric value indicating the line width used to plot the pbc box.}
  \item{cex.xyz}{a numeric value indicating the magnification used to plot the labels of the x, y and z axes.}
  \item{cex.abc}{a numeric value indicating the magnification used to plot the labels of the a, b and c axes.}
  \item{col}{a vector indicating the colors to use to plot each atom (If \code{length(col)!=length(xyz)/3} \code{col} is recycled). By default the Red Green Blue values corresponding to the atomic symbols are taken from the \code{elements} data set.}
  \item{radii}{either a character string indicating the type of radii or a numeric vector specifying the radii of each atom to use to plot atoms as spheres (see details).}
  \item{bg.col}{the color of the background of the scene.}
  \item{add}{a logical value indicating whether the plot has to be added to a existing scene (see \code{rgl.cur} and \code{open3d}).}
  \item{windowRect}{a vector of four integers indicating the left, top, right and bottom of the displayed window in pixels (see \code{par3d}).}
  \item{userMatrix}{a 4 by 4 matrix describing user actions to display the scene (see \code{par3d}).}
  \item{FOV}{the field of view. This controls the degree of parallax in the perspective view (see par3d).}
  \item{file}{the name of a PDB file to visualize.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
Three different visualization modes are allowed.
\itemize{
  \item   When \code{type="p"}: Points are drawn at each atomic positions (very light visualization mode).
  \item   When \code{type="l"}: Lines are drawn between bonded atoms. The connectivity of the system has to be specifyed.
  \item   When \code{type="s"}: Spheres are drawn at each atomic positions (heavy visualization mode).
  
  The radii of the spheres are given by \code{radii}.
  \itemize{
    \item   When \code{radii="rcov"}: Covalent radii, taken from the \code{elements} data set, are used.
    \item   When \code{radii="rvdw"}: Van der Waals radii, taken from the \code{elements} data set, are used (default).
    \item   When \code{radii} is a numeric vector: The numeric values are used to assign to each atom a radius. If \code{length(radii) != length(xyz)/3} \code{radii} is recycled.}
        }
Using \code{type="ls"} is allowed for "balls and sticks" representation. In that case, \code{radii} is set to half the covalent radii.
}
\value{
Return (using invisible) a two-column data.frame containing the IDs
and type indicators of the objects added to the scene.
}
\author{Julien Idé}

\seealso{
\code{\link{read.pdb}}, \code{\link{pymol}},\code{rgl.pop}, \code{rgl.cur}, \code{par3d}
}
\examples{
## Read a PDB file from the RCSB online database
pdb <- read.pdb("1bg2")
visualize(pdb, type = "p" , windowRect=c(0  ,0  ,400,400))
visualize(pdb, type = "l" , windowRect=c(400,0  ,800,400))
visualize(pdb, type = "s" , windowRect=c(0  ,400,400,800), radii = "rvdw")
visualize(pdb, type = "ls", windowRect=c(400,400,800,800))
}
\keyword{ dynamic }
