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
  plot.cna(x, weights=E(x$community.network)$weight*10,
           layout=NULL, col=NULL, ...)

  plot3d.cna(x, pdb = NULL, node.size = NULL,
             weights=E(x$community.network)$weight*10,
             layout = layout.pdb(pdb, x), col = NULL, ...)
}
\arguments{
  \item{x}{ A protein network graph object as obtained from the
    \sQuote{cna} function. }
  \item{community}{ A communities object as obtained from the
    \sQuote{cna} function and contained within x or a numeric vector of
    community or cluster membership. }
  \item{weights}{ A numeric vector containing the edge weights for the
    network. }
  \item{layout}{ Either a function or a numeric matrix. It specifies how the
    vertices will be placed on the plot. }
  \item{col}{ A vector of colors used for node/vertex rendering. If NULL
    these values are taken from the input network
    V(x$clustered.network)$color }
  \item{pdb}{ A PDB structure object obtained from \sQuote{read.pdb}. }
  \item{node.size}{ A numeric vector of node/community sizes. If NULL
    the size will be proportional to the residue membership of each community. }
  \item{\dots}{ Additional graphical parameters for
    \sQuote{plot.igraph}. }
}
\details{
  This function calls  \sQuote{plot.igraph} or \sQuote{plot.communities}
  from the igraph package to plot cna networks the way we like them.

  The plot layout is user setable, we like the options of: 
  \sQuote{layout.fruchterman.reingold},  \sQuote{layout.mds} or
  \sQuote{layout.svd}.

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

# Plot full residue network colored by communities - will be very slow!!
## plot.cna(net$network, net$communities) ## BUG!!!

# Play with plot settings
#plot.cna(net, layout=layout.mds(net$community.network), vertex.size=15,
#col=c("blue","green") ) ## BUG ignores vertex.color!!

\dontrun{
## 3D plot with RGL
plot3d.cna(net, pdb)
# rglplot(net$community.network)
}

}
\keyword{ hplot }
