\name{plot.dccm}
\alias{plot.dccm}
\title{ DCCM Plot }
\description{
  Plot a dynamical cross-correlation matrix.
}
\usage{
\method{plot}{dccm}(x, sse=NULL, colorkey=TRUE,
                    at=c(-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1),
                    main="Residue Cross Correlation",
                    helix.col = "gray20", sheet.col = "gray80",
                    inner.box=TRUE, outer.box=FALSE,
                    xlab="Residue No.", ylab="Residue No.",
                    margin.segments=NULL, segment.col=vmd.colors(), segment.min=1, ...)
}
\arguments{
  \item{x}{ a numeric matrix of atom-wise cross-correlations as output by
    the \sQuote{dccm} function. }
  \item{sse}{ secondary structure object as returned from
    \code{\link{dssp}}, \code{\link{stride}} or \code{\link{read.pdb}}. }
  \item{colorkey}{ logical, if TRUE a key is plotted. }
  \item{at}{ numeric vector specifying the levels to be colored. }
  \item{main}{ a main title for the plot. }
  \item{helix.col}{ The colors for rectangles representing alpha helices. }
  \item{sheet.col}{ The colors for rectangles representing beta strands. }
  \item{inner.box}{ logical, if TRUE an outer box is drawn. }
  \item{outer.box}{ logical, if TRUE an outer box is drawn. }
  \item{xlab}{ a label for the x axis. } 
  \item{ylab}{ a label for the y axis. }
  \item{margin.segments}{ a numeric vector of cluster membership as 
     obtained from cutree() or other community detection method. This 
     will be used for bottom and left margin annotation. }
  \item{segment.col}{ a vector of colors used for each cluster group in 
     margin.segments. }
  \item{segment.min}{ a single element numeric vector that will cause 
    margin.segments with a length below this value to be excluded from
    the plot. }
   \item{\dots}{ additional graphical parameters for contourplot. }
}
\details{
  See the \sQuote{contourplot} function from the lattice package for plot customization options, and the functions \code{\link{dssp}} and \code{\link{stride}} for further details.
}
\value{
  Called for its effect.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ Be sure to check the correspondence of your \sQuote{sse} object with the
  \sQuote{cij} values being plotted as no internal checks are currently performed.
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
  inds <- atom.select(pdb, resno=c(24:27,85:90), elety='CA')

  ## lsq fit of trj on pdb
  xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

  ## Dynamic cross-correlations of atomic displacements
  cij <- dccm(xyz)

  ## Default plot
  plot.dccm(cij)

  ## Change the color scheme and the range of colored data levels
  plot.dccm(cij, contour=F, col.regions=bwr.colors(200), at=seq(-1,1,by=0.01) )

  ## Add secondary structure annotation to plot margins
  sse <- dssp(read.pdb("1W5Y"), resno=FALSE)
  plot.dccm(cij, sse=sse) 

  ## Add additional margin annotation for chains..
  ch <- ifelse(pdb$atom[pdb$calpha,"chain"]=="A", 1,2)
  plot.dccm(cij, sse=sse, margin.segments=ch)

  ## Plot with cluster annotation from dynamic network analysis
  #net <- cna(cij)
  #plot.dccm(cij, margin.segments=net$raw.communities$membership)

  ## Focus on major communities (i.e. exclude those below a certain total length)
  #plot.dccm(cij, margin.segments=net$raw.communities$membership, segment.min=25)

}

}
\keyword{ hplot }
