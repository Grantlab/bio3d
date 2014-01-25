\name{identify.cna}
\alias{identify.cna}
\title{
 Label Idetification in CNA Plots
}
\description{
  This function add labels to the objects in a two-dimensional CNA plot.
}
\usage{
  \method{identify}{cna}(x, labels=NULL, cna=NULL, labelshorten=TRUE, ...)
}
\arguments{
  \item{x}{ A numeric matrix with Nx2 dimensions, where N is equal to
    the number of objects in a 2D CNA plot such as obtained from
    \sQuote{plot.cna} function. }
  \item{labels}{ An optional character vector giving labels for the points.
    Will be coerced using \sQuote{as.character}, and recycled if
    necessary to the length of \sQuote{x}.  Excess labels will be
    discarded, with a warning. }
  \item{cna}{ A network object as returned from the \sQuote{cna} function. }
  \item{labelshorten}{ Logical, if TRUE an attempt will be made to shorten
    the printed and returned lables. }
  \item{\dots}{ Extra options passed to \sQuote{identify} function. } 
}
\details{
  This function calls the \sQuote{identify} function and adds labels to the objects
  in a two-dimensional CNA plot by left clicking over them with the
  mouse. A click with the right mouse button will stop the function.
}
\value{
  If \sQuote{labels} or \sQuote{cna} inputs are provided then a membership vector 
  will be returned with the selected communitiy ids and their members. Otherwise a 
  vector with the ids of the selected communities will be returned.
}
\author{ Guido Scarabelli and Barry Grant }

\seealso{ \code{\link{identify}}, \code{\link{plot.igraph}}, \code{\link{plot.communities}}, \code{\link{igraph.plotting}} \code{\link{plot.cna}}}
\examples{

# Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb,"///24:27,85:90///CA/")

# lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

# calculate dynamical cross-correlation matrix
cij <- dccm(xyz)

# Build, and betweeness cluster, a network graph 
net <- cna(cij)

# Plot the network
xy <- plot.cna(net)

# Use identify.cna on the communities
d <- identify.cna(xy, cna=net)

# Right click to end the function...
  ## d <- identify.cna(xy, labels=net$communities$membership)
  ## d <- identify(xy, cnet=net)
  ## d <- identify(xy, summary(net)$members)

}
\keyword{ utility }