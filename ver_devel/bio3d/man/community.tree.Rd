\name{community.tree}
\alias{community.tree}
\title{
 Reconstruction of the Girvan-Newman community tree for a cna class object.
}
\description{
  This function reconstructs the community tree of the community
  analysis performed by \sQuote{cna} function. It allows the user to
  explore different community partitionigs for the network.
}
\usage{
  community.tree(x, rescale=FALSE)
}
\arguments{
  \item{x}{  A protein network graph object as obtained from the
    \sQuote{cna} function. }
  \item{rescale}{ Logical, indicating whether to rescale the community
    names starting from 1. If FALSE, the community names will start from
    N+1, where N is the number of nodes. } 
}
\value{
  Returns a list object that includes the following components:
  \item{modularity}{ A numeric vector containing the modularity values. }
  \item{tree}{ A numeric matrix containing in each row the
    community residue memberships corresponding to a modularity
    value. The rows are ordered according to the \sQuote{modularity} object. }
  \item{num.of.comms}{ A numeric vector containing the number of
    communities per modularity value. The vector elements are ordered accoring
    to the \sQuote{modularity} object. }
}
\details{
  The input of this function should be a \sQuote{cna} class object
  containing a $network, a $communities$merges and a $modularity
  attributes. 
  The function reconstructs the community residue memberships for each
  modularity value. The purpose is to identify the most significant
  community partitioning, which it might be obtained with a value close
  to the maximum of the modularity, but not with the maximum itself.
}  
\seealso{
   \code{\link{cna}}, \code{\link{summary.cna}}
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

tree.net <- community.tree(net, rescale=TRUE)

# Get the maximum modularity
max.modularity <- which(tree.net$modularity ==
max(tree.net$modularity))

# Same as the communities in net object
tree.net$tree[max.modularity,]

# New membership partitioning
tree.net$tree[max.modularity+2,]

}
\keyword{analysis}