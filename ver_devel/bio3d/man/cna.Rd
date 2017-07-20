\name{cna}
\alias{cna}
\alias{cna.dccm}
\alias{cna.ensmb}
\title{ Protein Dynamic Correlation Network Construction and Community Analysis.  }
\description{
  This function builds both residue-based and community-based undirected weighted
  network graphs from an input correlation matrix, as obtained from the
  functions \sQuote{dccm}, \sQuote{dccm.nma}, and \sQuote{dccm.enma}. Community
  detection/clustering is performed on the initial residue based 
  network to determine the community organization and network structure
  of the community based network.
}
\usage{
  cna(cij, \dots)
  \method{cna}{dccm}(cij, cutoff.cij=0.4, cm=NULL,  vnames=colnames(cij),
      cluster.method="btwn", collapse.method="max", 
      cols=vmd_colors(), minus.log=TRUE, \dots)
  \method{cna}{ensmb}(cij, \dots, ncore = NULL)
}
\arguments{
  \item{cij}{ A numeric array with 2 dimensions (nXn) containing atomic
    correlation values, where "n" is the residue number. The matrix
    elements should be in between 0 and 1 (atomic correlations). Can be
    also a set of correlation matrices for ensemble network analysis. See
     \sQuote{dccm} function in bio3d package for further details. }
  \item{\dots}{  Additional arguments passed to the methods
    \code{cna.dccm} and \code{cna.ensmb}. }
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
    btwn="Girvan-Newman betweenness"\cr
    walk="Random walk"\cr
    greed="Greedy algorithm for modularity optimization"\cr 
    infomap="Infomap algorithm for community detection"\cr }  
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
  \item{ncore}{ Number of CPU cores used to do the calculation. By default,
    use all available cores. }
}
\value{
  Returns a list object that includes igraph network and community
  objects with the following components:
  \item{network}{ An igraph residue-wise graph object. See below
    for more details.}
  \item{communities}{ An igraph residue-wise community object. See
    below for more details. }
  \item{communitiy.network}{ An igraph community-wise graph object. See
    below for more details. }
  \item{community.cij}{ Numeric square matrix containing the absolute values of
    the atomic correlation input matrix for each community as obtained
    from \sQuote{cij} via application of  \sQuote{collapse.method}. }
  \item{cij}{ Numeric square matrix containing the absolute values of
    the atomic correlation input matrix. }
}
\details{
  The input to this function should be a correlation matrix as
  obtained from the  \sQuote{dccm}, \sQuote{dccm.mean} or \sQuote{dccm.nma}
  and related functions. Optionally, a contact map \sQuote{cm} may also given 
  as input to filter the correlation matrix resulting in the exclusion of  
  network edges between non-contacting atom pairs (as defined in the contact map).

  Internally this function calls the igraph package
  functions \sQuote{graph.adjacency},
  \sQuote{edge.betweenness.community}, \sQuote{walktrap.community},
  \sQuote{fastgreedy.community}, and \sQuote{infomap.community}. The first constructs an undirected
  weighted network graph. The second performs Girvan-Newman style clustering by
  calculating the edge betweenness of the graph, removing the edge with
  the highest edge betweenness score, calculates modularity (i.e. the difference between
  the current graph partition and the partition of a random graph, see Newman and Girvan,
  Physical Review E (2004), Vol 69, 026113), then recalculating edge betweenness
  of the edges and again removing the one with the highest score, etc. 
  The returned community partition is the one with the highest overall modularity value.
  \sQuote{walktrap.community} implements the Pons and Latapy algorithm based on the idea that
  random walks on a graph tend to get "trapped" into densely connected parts of it,
  i.e. a community. The random walk process is used to determine a distance
  between nodes. Nodes with low distance values are joined in the same community.
  \sQuote{fastgreedy.community} instead determines the community structure based on the
  optimization of the modularity. In the starting state each node is isolated and 
  belongs to a separated community. Communities are then joined together (according to the
  network edges) in pairs and the modularity is calculated. At each step the join resulting
  in the highest increase of modularity is chosen. This process is repeated until a single
  community is obtained, then the partitioning with the highest modularity score is selected.
  \sQuote{infomap.community} finds community structure that minimizes the expected description
     length of a random walker trajectory.
}
\author{ Guido Scarabelli and Barry Grant }  
\seealso{
   \code{\link{plot.cna}},  \code{\link{summary.cna}},
   \code{\link{vmd.cna}}, \code{\link[igraph:graph.adjacency]{graph.adjacency}},
   \code{\link[igraph:edge.betweenness.community]{edge.betweenness.community}},
   \code{\link[igraph:walktrap.community]{walktrap.community}}, 
   \code{\link[igraph:fastgreedy.community]{fastgreedy.community}},
   \code{\link[igraph:infomap.community]{infomap.community}}
}
\examples{
\donttest{
# PDB server connection required - testing excluded

if (!requireNamespace("igraph", quietly = TRUE)) {
   message('Need igraph installed to run this example')
} else {

##-- Build a correlation network from NMA results
## Read example PDB
pdb <- read.pdb("4Q21")

## Perform NMA
modes <- nma(pdb)
#plot(modes, sse=pdb)

## Calculate correlations 
cij <- dccm(modes)
#plot(cij, sse=pdb)

## Build, and betweenness cluster, a network graph
net <- cna(cij, cutoff.cij=0.35)
#plot(net, pdb)

## within VMD set 'coloring method' to 'Chain' and 'Drawing method' to Tube
#vmd.cna(net, trim.pdb(pdb, atom.select(pdb,"calpha")), launch=TRUE )

##-- Build a correlation network from MD results
## Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb, resno=c(24:27,85:90), elety='CA')

## lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

## calculate dynamical cross-correlation matrix
cij <- dccm(xyz)

## Build, and betweenness cluster, a network graph
net <- cna(cij)

# Plot coarse grained network based on dynamically coupled communities
xy <- plot.cna(net)
plot.dccm(cij, margin.segments=net$communities$membership)


##-- Begin to examine network structure - see CNA vignette for more details
net
summary(net)
attributes(net)
table( net$communities$members )

}
}
}
\keyword{analysis}
