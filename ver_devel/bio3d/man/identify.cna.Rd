\name{identify.cna}
\alias{identify.cna}
\title{
 Identify Points in a CNA Protein Structure Network Plot
}
\description{
  \sQuote{identify.cna} reads the position of the graphics pointer when the
  (first) mouse button is pressed.  It then searches the coordinates given
  in \sQuote{x} for the point closest to the pointer.  If this point is 
  close enough to the pointer, its index and community members will be
  returned as part of the value of the call and the community members will 
  be added as labels to the plot.
}
\usage{
  \method{identify}{cna}(x, labels=NULL, cna=NULL, ...)
}
\arguments{
  \item{x}{ A numeric matrix with Nx2 dimensions, where N is equal to
    the number of objects in a 2D CNA plot such as obtained from the 
    \sQuote{plot.cna} and various \sQuote{layout} functions. }
  \item{labels}{ An optional character vector giving labels for the points.
    Will be coerced using \sQuote{as.character}, and recycled if
    necessary to the length of \sQuote{x}.  Excess labels will be
    discarded, with a warning. }
  \item{cna}{ A network object as returned from the \sQuote{cna} function. }
  \item{\dots}{ Extra options passed to \sQuote{identify} function. } 
}
\details{
  This function calls the \sQuote{identify} and \sQuote{summary.cna} functions to 
  query and label 2D CNA protein structure network plots produced by the \sQuote{plot.cna} function. Clicking with the mouse on plot points will add the corresponding  
  labels and them to the plot and returned list object. A click with the right 
  mouse button will stop the function.
}
\value{
  If \sQuote{labels} or \sQuote{cna} inputs are provided then a membership vector 
  will be returned with the selected community ids and their members. Otherwise a 
  vector with the ids of the selected communities will be returned.
}
\author{ Guido Scarabelli and Barry Grant }

\seealso{ \code{\link{plot.cna}}, \code{\link{identify}}, \code{\link{plot.igraph}}, \code{\link{plot.communities}}, \code{\link{igraph.plotting}} }
\examples{
\dontrun{
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
## d <- identify(xy, summary(net)$members)
}
}
\keyword{ utility }