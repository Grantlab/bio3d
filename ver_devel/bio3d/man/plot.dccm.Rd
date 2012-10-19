\name{plot.dccm}
\alias{plot.dccm}
\title{ DCCM Plot }
\description{
  Plot a dynamical cross-correlation matrix.
}
\usage{
plot.dccm(x, sse = NULL, colorkey = TRUE, at = c(-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1), main = "Residue Cross Correlation", helix.col = "gray20", sheet.col = "gray80", inner.box = TRUE, outer.box = FALSE, xlab = "Residue No.", ylab = "Residue No.", ...)
}
\arguments{
  \item{x}{ a numeric matrix of atom-wise cross-correlations as output by
    the \sQuote{dccm} function. }
  \item{sse}{ secondary structure object as returned from
    \code{\link{dssp}} or \code{\link{stride}}. }
  \item{colorkey}{ logical, if TRUE a key is plotted. }
  \item{at}{ numeric vector specifying the levels to be colored. }
  \item{main}{ a main title for the plot. }
  \item{helix.col}{ The colors for rectangles representing alpha helices. }
  \item{sheet.col}{ The colors for rectangles representing beta strands. }
  \item{inner.box}{ logical, if TRUE an outer box is drawn. }
  \item{outer.box}{ logical, if TRUE an outer box is drawn. }
  \item{xlab}{ a label for the x axis. } 
  \item{ylab}{ a label for the y axis. }
  \item{\dots}{ additional graphical parameters for image. }
}
\details{
  See the functions \sQuote{plot.default}, \code{\link{dssp}} and
  \code{\link{stride}} for further details.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ Be sure to check the correspondence of your \sQuote{sse} object with the
  \sQuote{cij} values being plotted as no internal checks are
  performed.
}
\seealso{ \code{\link{plot.bio3d}}, \code{\link{plot.dmat}},
  \code{\link{filled.contour}}, \code{\link{contour}},
  \code{\link{image}} \code{\link{plot.default}}, \code{\link{dssp}},
  \code{\link{stride}} }
\examples{

\dontrun{
  ##-- Read example trajectory file
  trtfile <- system.file("examples/hivp.dcd", package="bio3d")
  trj <- read.dcd(trtfile)

  ## Read the starting PDB file to determine atom correspondence
  pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
  pdb <- read.pdb(pdbfile)

  ## select residues 24 to 27 and 85 to 90 in both chains
  inds <- atom.select(pdb,"///24:27,85:90///CA/")

  ## lsq fit of trj on pdb
  xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

  cij <- dccm(xyz)
  plot.dccm(cij)

  ## SSE axis labs
  sse <- dssp(pdb, resno=FALSE)
  gridpoints <- c( sse$helix$start[sse$helix$length>4],
                  sse$sheet$start[sse$sheet$length>4])
  
  scales <- list(at=gridpoints,
                 labels=gridpoints, alternating=1, ticks=F)

  plot.dccm(cij, sse=sse, scales=scales)
  ##add.grid(sse)

  ## Plot with changed colors
  plot.dccm(cij, col.regions=rev(heat.colors(12)))

}

}
\keyword{ hplot }
