\name{network.amendment}
\alias{network.amendment}
\title{
  Amendment of cna Class Object After Provinding a New Community Membership Vector.
}
\description{
  This function changes the attributes $communities$membership,
  $community.network and $community.cij for a \sQuote{cna} class object.
}
\usage{
  network.amendment(x, membership, minus.log=TRUE)

}
\arguments{
  \item{x}{  A protein network graph object as obtained from the
    \sQuote{cna} function. }
  \item{membership}{ A numeric vector containing the new community
    membership. } 
  \item{minus.log}{ Logical. Whether to use the minus.log on the cij
    values. }
}
\value{
  Returns a \sQuote{cna} class object with the attributes changed
  according to the membership vector provided.
}
\details{
  This function is useful to change the community partitioning of a
  \sQuote{cna} object, according to the community membership vector
  provided (i.e. it should be used with the output of the \sQuote{community.tree} function).
}  
\seealso{
   \code{\link{cna}}, \code{\link{summary.cna}}, \code{\link{community.tree}}
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
tree.net$tree[max.modularity+3,]

# New cna object
net.optimized <- network.amendment(net, tree.net$tree[max.modularity+3,], minus.log=TRUE)

par(mfcol=c(1,2))
plot.cna(net)
plot.cna(net.optimized)

#view.cna(net.optimized, pdb, launch=TRUE)
}
\keyword{utility}