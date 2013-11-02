\name{cna-package}
\alias{cna-package}
\alias{cna}
\docType{package}
\title{
Correlation Network Analysis
}
\description{
This package has utilities to carry out basic analysis on the correlation networks.
It is integrated with bio3d, igraph, graph R packages. It also reads the output of the
VMD-DNA plugin analyisis, and compare them with the ones from the network R packages.
Below you can find the function list.

Cna contains the following functions:
- dccm.mean (average 3D cij.matri\cr
- cmap.filter (average 3D contact.matrix)\cr
- igraph.btwn (betweenness values with igraph)\cr
- vmd.btwn (betweenness values of vmd-plugin)\cr
- cna (community clustering with igraph)\cr
- plot.cna (plot cna graphs)\cr
- layout.pdb (generates x,y,z coordinates to plot communities)\cr
- identify.cna (equivalent of \sQuote{identify} function for cna plots)\cr
- plot3d.cna (plot cna communities in a 3D plot)\cr
- print.cna (print cna function output)\cr
- prune.cna (delete specified communities)\cr
- critical residues (critical nodes list)\cr
- membership.comparison (plot the community membership of specified residues)\cr
- vmd.comms (vmd-plugin communities and critical nodes)\cr
- cluster.sim (jaccard, rand and adjusted rand index calculation)\cr
- dssp.block (coarse grain the cij plot according to the sse)\cr

Below there are three examples, that should guide the user through the analysis available in the package.
The examples work in a consecutive way, so you need inputs in example 1 to run example 2, and so on...
All the examples shown use a simple 8 node network, but other three weighted network (with 10, 50 and 100 nodes)
are available with "load(paste(system.file("data",package="cna"),"/example.node10.RData",sep=""))",
"data(node.50)" or "data(node.100)".\cr
For VMD-DNA-plugin, the files "vmd-50-node.out", "vmd-btwn-50-node.dat",
"vmd-100-node.out", "vmd-btwn-100-node.dat" are provided in "cna/data" folder.\cr
Upon request, the pdb and trajectory data from 10 kinesin APO cMD simulations (kinesin.pdb, kin.trj1, kin.trj2, kin.trj3, etc...)
are available (not included due to the MBs required).

}
\details{
\tabular{ll}{
Package: \tab cna\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-12\cr
License: \tab What license is it under?\cr
}
}

\keyword{package}
\seealso{
dccm (bio3d), cmap (bio3d), 
graph.adjacency (igraph), edge.betweenness (igraph), edge.betweenness.community (igraph), 
graphAM (graph), brandes.betweenness.centrality (graph)
}


\examples{
\dontrun{
## Example 1 -- Network generation
  ## igraph package
  library(igraph)

  ## adjacency matrix for 8 nodes
  data(node.8) 

  ## Network geneation, inputs: adjacency matrix (AM), edge type (directed/undirected), weights (taken from the AM), values for the diagonal.
  network <- graph.adjacency(node.8,mode="undirected",weighted=TRUE,diag=FALSE)	

  ## prints out the number of nodes and edges
  network    

  ## list of the nodes
  V(network) 

  ## list of the edges
  E(network) 

  ## generates a 2D plot of the network
  plot(network)  



  ## graph package
  library(graph)
  
  ## adjacency matrix for 8 nodes
  data(node.8) 

  ## Network generation, inputs: adjacency matrix (AM), edge type (directed/undirected), edge weights defalut values (needed to use AM(ij) as weights!)
  network <- graphAM(adjMat=node.8, edgemode="undirected",values=list(weight=1))	
 
  ## prints out the number of nodes and edges
  network

  ## prints out the edge weights
  edgeWeights(network)

  ## generates a 2D plot of the network
  plot(network)  

######################################

## Example 2 -- Community calculation
  ## igraph package -- igraph community object
  net <- cna(node.8)

  ## igraph community membership
  net$raw.network$membership
  plot.cna(net)
  
  ## vmd community list
  vmd.comms.rslt <- vmd.comms(paste(pth,"/vmd-8-node.out",sep=""),adjmatrix=node.8)
  plot(vmd.comms.rslt$network,edge.witdh=vmd.comms.rslt$weights)
  

######################################

## Example 3 -- Check community consistency
  comparison <- compare.communities(net$raw.communities,vmd.comms.rslt$membership,method="rand")

######################################

}
}
