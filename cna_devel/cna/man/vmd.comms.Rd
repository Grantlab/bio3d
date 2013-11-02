\name{vmd.comms}
\alias{vmd.comms}
\title{
VMD Community Structure
}
\description{
This function reads a VMD-DNA-plugin community output files and creates a numeric residue community-membership vector.
}
\usage{
vmd.comms(vmd.out.file=vmd.out.file,adjmatrix=adjmatrix,vnames=colnames(adjmatrix))
}
\arguments{
  \item{vmd.out.file}{
Community output file from VMD-DNA-plugin. This file is generated from a weighted adjacency matrix obtained with the multiplication of
an atomic correlation matrix (cij matrix) with a contact matrix (cm matrix). The gncommunity script (implemented with the girvan-newman
algorithm) is used to carry out the analysis.
}
  \item{adjmatrix}{
A numeric square adjacency matrix object, as obtained from adjmatrix.setup.
}
  \item{vnames}{
A vector of names for each column in the input adjmatrix. 
}  
}
\value{
This function reads the output of VMD-DNA analysis and make it compatible with igraph object format.
An object containing the following attributes is created:\cr
- $network (an igraph network object)\cr
- $membership (a numeric vector of length n (where "n" is the number of residues) with the community-membership)\cr
- $weights (a numeric vector with the edge weights for the plot)\cr
- $critical.nodes (a numeric table containing the information for each critical. Each row show the betweenness value of that edge (column 1),
  the nodes belonging to that edge (columns 2 and 3) and their community community memberships (columns 4 and 5).
}

\examples{
## path for the input
pth <- system.file("data",package="cna")  ## path to the example files

## load adjmatrix
data(node.8)

## vmd residue community-membership vector
vmd.comms.rslt <- vmd.comms(paste(pth,"/vmd-8-node.out",sep=""),adjmatrix=node.8)

## calling igraph package for plotting
library(igraph)
plot(vmd.comms.rslt$network,  edge.width=vmd.comms.rslt$weights)

}
\keyword{utility}

