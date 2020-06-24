\name{plot.cna}
\alias{plot.cna}
\title{ Protein Structure Network Plots in 2D and 3D.  }
\description{
  Plot a protein dynamic network as obtained from the \emph{cna} function.
}
\usage{
  \method{plot}{cna}(x, pdb = NULL, weights=NULL, vertex.size=NULL,
            layout=NULL, col=NULL, full=FALSE, scale=TRUE, color.edge = FALSE, 
            interactive=FALSE, ...)
}
\arguments{
  \item{x}{ A protein network graph object as obtained from the
    \sQuote{cna} function. }
  \item{pdb}{ A PDB structure object obtained from \sQuote{read.pdb}. 
    If supplied this will be used to guide the network plot \sQuote{layout}, 
    see \sQuote{layout.cna} for details. }    
  \item{weights}{ A numeric vector containing the edge weights for the
    network. }
  \item{vertex.size}{ A numeric vector of node/community sizes. If NULL
    the size will be taken from the input network graph object \sQuote{x}. 
    Typically for \sQuote{full=TRUE} nodes will be of an equal size and 
    for \sQuote{full=FALSE} community node size will be proportional to 
    the residue membership of each community. }
  \item{layout}{ Either a function or a numeric matrix. It specifies how the
    vertices will be placed on the plot. See \sQuote{layout.cna}. }
  \item{col}{ A vector of colors used for node/vertex rendering. If NULL
    these values are taken from the input network
    \sQuote{V(x$community.network)$color}. }
  \item{full}{ Logical, if TRUE the full all-atom network rather than the 
    clustered community network will be plotted. }    
  \item{scale}{ Logical, if TRUE weights are scaled with respect to the network. }   
  \item{color.edge}{ Logical, if TRUE edges are colored with respect to their weights. }   
  \item{interactive}{ Logical, if TRUE interactive graph will be drawn where users
    can manually adjust the network (positions of vertices, colors of edges, etc.).
    Needs Tcl/Tk support in the installed R build. }
  \item{\dots}{ Additional graphical parameters for \sQuote{plot.igraph}. }
}
\details{
  This function calls \sQuote{plot.igraph} from the igraph package to 
  plot cna networks the way we like them.

  The plot layout is user settable, we like the options of: 
  \sQuote{layout.cna}, \sQuote{layout.fruchterman.reingold},  
  \sQuote{layout.mds} or \sQuote{layout.svd}. Note that first of these 
  uses PDB structure information to produce a more meaningful layout.

  Extensive plot modifications are possible by setting additional graphical 
  parameters (\dots). These options are detailed in \sQuote{igraph.plotting}. 
  Common parameters to alter include:

  \describe{
    \item{vertex.label:}{Node labels, \code{V(x$network)$name}. Use NA to omit.}
    \item{vertex.label.color:}{Node label colors, see also \code{vertex.label.cex} etc. }
    \item{edge.color:}{Edge colors,  \code{E(x$network)$color}. }
    \item{mark.groups:}{Community highlighting, a community list object, see also \code{mark.col} etc.}
  }
  
}
\value{
  Produces a network plot on the active graphics device.
  Also returns the plot layout coordinates silently, which can be passed to the 
  \sQuote{identify.cna} function.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant and Guido Scarabelli }
\note{ Be sure to check the correspondence of your \sQuote{pdb} object with your
  network object \sQuote{x}, as few internal checks are currently performed by the 
  \sQuote{layout.cna} function.
}
\seealso{ \code{\link[igraph:plot.igraph]{plot.igraph}}, 
  \code{\link[igraph:communities]{plot.communities}},
  \code{\link[igraph:plot.common]{igraph.plotting}} }
\examples{
\donttest{
# PDB server connection required - testing excluded

if (!requireNamespace("igraph", quietly = TRUE)) {
   message('Need igraph installed to run this example')
} else {

##-- Build a CNA object
pdb <- read.pdb("4Q21")
modes <- nma(pdb)
cij <- dccm(modes)
net <- cna(cij, cutoff.cij=0.2)

# Plot coarse grain network based on dynamically coupled communities
xy <- plot.cna(net)
#plot.dccm(cij, margin.segments=net$communities$membership)

# Chose a different PDB informed layout for plot
plot.cna(net, pdb)

# Play with plot layout and colors...
plot.cna(net, layout=igraph::layout.mds(net$community.network), col=c("blue","green") ) 

# Plot full residue network colored by communities - will be slow due to number of edges!!
plot.cna(net, pdb, full=TRUE)

# Alter plot settings
plot.cna(net, pdb, full=TRUE, vertex.size=3, weights=1, vertex.label=NA)

}
}
}
\keyword{ hplot }
