\name{plot.cna}
\alias{plot.cna}
\alias{plot3d.cna}
\title{
 Protein Dynamic Network Plots in 2D and 3D.
}
\description{
  Plot a protein dynamic network as obtained from the \sQuote{cna} function.
}
\usage{
  \method{plot}{cna}(x, weights=NULL,
           layout=NULL, col=NULL, full=FALSE, ...)

  plot3d.cna(x, pdb = NULL, node.size = NULL,
             weights=NULL,
             layout = layout.cna(x, pdb, k=3), col = NULL, ...)
}
\arguments{
  \item{x}{ A protein network graph object as obtained from the
    \sQuote{cna} function. }
  \item{weights}{ A numeric vector containing the edge weights for the
    network. }
  \item{layout}{ Either a function or a numeric matrix. It specifies how the
    vertices will be placed on the plot. See \sQuote{layout.cna}. }
  \item{col}{ A vector of colors used for node/vertex rendering. If NULL
    these values are taken from the input network
    V(x$community.network)$color }
  \item{full}{ Logical, if TRUE the full all-atom network rather than the 
    clustered community  network will be plotted. }       
  \item{pdb}{ A PDB structure object obtained from \sQuote{read.pdb}. }
  \item{node.size}{ A numeric vector of node/community sizes. If NULL
    the size will be proportional to the residue membership of each community. }
  \item{\dots}{ Additional graphical parameters for \sQuote{plot.igraph}. }
}
\details{
  This function calls \sQuote{plot.igraph} or \sQuote{plot.communities}
  from the igraph package to plot cna networks the way we like them.

  The plot layout is user settable, we like the options of: 
  \sQuote{layout.cna}, \sQuote{layout.fruchterman.reingold},  
  \sQuote{layout.mds} or \sQuote{layout.svd}. Note that first of these 
  uses PDB structure information to produce a more meaningful layout.

  \dots can contain:
  col=vmd.colors(),
  mark.col=vmd.colors(alpha=0.3),
  mark.border=vmd.colors()
  vertex.size:  Node sizes:   V(x$network)size
  vertex.color: Node colors:  V(x$network)color
  vertex.label: Node labels:  vertex ids
  edge.width:   Edge weights: E(x$network)$weight
  edge.color:   Edge colors:  E(x$network)$color
      (also vertex.label.color, vertex.label.cex etc.
       see \code{\link{igraph.plotting}} for others.)
  
}
\value{
  Called for its effect.
  Also returns the plot layout coordinates silently.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{ Be sure to check the correspondence of your \sQuote{sse} object with the
  \sQuote{cij} values being plotted as no internal checks are
  performed.
}
\seealso{ \code{\link{plot.igraph}}, \code{\link{plot.communities}},
  \code{\link{igraph.plotting}} }
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
#plot(cij)

# Build, and betweeness cluster, a network graph 
net <- cna(cij)

# Plot coarse grain network based on dynamicaly coupled communities
xy <- plot.cna(net)
plot.dccm2(cij, margin.segments=net$communities$membership)

# Plot the community network
plot.cna(net)

# Chose a different PDB informed layout
plot.cna(net, layout=layout.cna(net,pdb))

# Plot full residue network colored by communities - will be very slow!!
plot.cna(net, layout=layout.cna(net,pdb, full=T), full=T, vertex.size=10, weights=1)


# Play with plot settings example
plot.cna(net, layout=layout.mds(net$community.network), col=c("blue","green") ) 

\dontrun{
## 3D plot with RGL
plot3d.cna(net, pdb)
# rglplot(net$community.network)
}

}
\keyword{ hplot }
