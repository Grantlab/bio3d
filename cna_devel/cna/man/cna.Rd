\name{cna}
\alias{cna}
\title{
  Protein Dynamic Correlation Network Construction and Community Analysis
}
\description{
  This function builds both residue-based and community-based undirected weighted
  network graphs from an input correlation matrix, as obtained from the
  functions \sQuote{dccm} and \sQuote{dccm.nma}. Community
  detection/clustering is performed on the initial residue based (raw)
  network to determine the community organization and network structure
  of the community based (clustered) network.
}
\usage{
  cna(cij, vnames=colnames(cij), cutoff.cij=0.4, cm=NULL, method="btwn")
  cna(cij, cutoff.cij=0.4, cm=NULL,  vnames=colnames(cij),
      cluster.method="btwn", collapse.method="max", 
      cols=vmd.colors(), minus.log=TRUE)
}
\arguments{
  \item{cij}{ A numeric array with 2 dimensions (nXn) containing atomic
    correlation values, where "n" is the residue number. The matrix
    elements should be in between 0 and 1 (atomic correlations). See
     \sQuote{dccm} function in bio3d package for further details. }
  \item{cutoff.cij}{ Numeric element specifying the cutoff on cij matrix
    values. Coupling below cutoff.cij are set to 0. }
  \item{cm}{ (optinal) A numeric array with 2 dimensions (nXn)
    containing binary contact values, where "n" is the residue
    number. The matrix elements should be 1 if two residues are in
    contact and 0 if not in contact. See the  \sQuote{cmap} function in
    bio3d package for further details. }
  \item{vnames}{ A vector of names for each column in the input
    cij. This will be used for referencing residues in a similar way
    to residue numbers in later analysis. }  
  \item{cluster.method}{ A character string specifying the method for
    community determination. Supported methods are:\cr
    btwn='Girvan-Newman betweenness'\cr
    walk='Random walk'\cr
    greedy='Greedy algorithm for modularity optimization'\cr }  
  \item{collapse.method}{ A single element character vector specifing
    the \sQuote{cij} collapse method, can be one of \sQuote{max}, \sQuote{median},
    \sQuote{mean}, or \sQuote{trimmed}. By defualt the \sQuote{max} method is
    used to collapse the input residue based \sQuote{cij} matrix into a smaller
    community based network by taking the maximium \sQuote{abs(cij)} value between
    communities as the comunity-to-community cij value for clustered
    network construction. }
  \item{cols}{ A vector of colors assigned to network nodes. }
  \item{minus.log}{ Logical, indicating whether \sQuote{-log(abs(cij))}
    values should be used for network construction. } 
}
\value{
  Returns a list object that includes igraph network and community
  objects with the following components:
  \item{raw.network}{  An igraph residue-wise graph object. See below
    for more details.}
  \item{raw.communities}{ An igraph residue-wise community object. See
    below for more details. }
  \item{clustered.network}{ An igraph community-wise graph object. See
    below for more details. }
  \item{clustered.communities}{ An igraph community-wise community
    object. See below for more details. }
  \item{clustered.cij}{ Numeric square matrix containing the absolute values of
    the atomic correlation input matrix for each community as obtained
    from \sQuote{raw.cij} via application of  \sQuote{collapse.method}. }
  \item{raw.cij}{ Numeric square matrix containing the absolute values of
    the atomic correlation input matrix. }
}
\details{
  The input to this function should be a correlation matrix as
  obtained from the  \sQuote{dccm}, \sQuote{dccm.mean} or \sQuote{dccm.nma}
  functions. Optionaly, a contact map may also given as input to
  filter the correlation matrix.

  TO BE UPDATED
  
  Internally this function calls the igraph package
  functions \sQuote{graph.adjacency},
  \sQuote{edge.betweenness.community}, \sQuote{walktrap.community},
  \sQuote{fastgreedy.community}. The first constructs an undirected
  weighted network graph. The second performs Girvan-Newman style clustering by
  calculating the edge betweenness of the graph, removing the edge with
  the highest edge betweenness score, calculates modularity (i.e. the difference between
  the current graph partition and the partition of a random graph, see Newman and Girvan,
  Physical Review E (2004), Vol 69, 026113), then recalculating edge betweenness
  of the edges and again removing the one with the highest score, etc. 
  The best community partition is the one with the highest modularity value.
  \sQuote{walktrap.community} implements the Pons and Latapy algorithm based on the idea that
  random walks on a graph tend to get “trapped” into densely connected parts of it,
  i.e. a community. The random walk process is used to determine a distance
  between nodes. Nodes with low distance values are joined in the same community.
  \sQuote{fastgreedy.community} instead determines the community structure based on the
  optimization of the modularity. In the starting state each node is isolated and 
  belongs to a separated community. Communities are then joined together (according to the
  network edges) in pairs and the modularity is calculated. At each step the join resulting
  in the highest increase of modularity is chosen. This process is repeated until a single
  community is obtained, then the partitioning with the highest modularity score is selected.
}  
\seealso{
   \code{\link{plot.cna}},  \code{\link{summary.cna}},
   \code{\link{view.cna}}, \code{\link{graph.adjacency}},
   \code{\link{edge.betweenness.community}},
  \code{\link{walktrap.community}}, \code{\link{fastgreedy.community}},
}
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
plot.dccm2(cij, margin.segments=net$raw.communities$membership)

# Plot full residue network colored by communities - will be very slow!!
## plot.cna(net$raw.network, net$raw.communities) ## BUG!!!

# Play with plot settings
#plot.cna(net, layout=layout.mds(net$clustered.network), vertex.size=15,
#col=c("blue","green") ) ## BUG ignores vertex.color!!

\dontrun{
## Example adjacency matrix input
data(node.50)
range(node.50)

comms <-cna(node.50,method="btwn",cutoff.cij=0.2)

# Plot colored by community membership
par(mfcol=c(1,2), mar=c(0,0,0,0))
plot(comms$raw.network)
plot(comms)

# Examine community structure
comms

attributes(comms)
comms$raw.communities$members
## More examples
}
}
\keyword{analysis}
