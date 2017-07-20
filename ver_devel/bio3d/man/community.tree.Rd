\name{community.tree}
\alias{community.tree}
\title{ Reconstruction of the Girvan-Newman Community Tree for a CNA Class Object.  }
\description{
  This function reconstructs the community tree of the community
  clustering analysis performed by the \sQuote{cna} function. 
  It allows the user to explore different network community partitions.
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
    communities per modularity value. The vector elements are ordered according
    to the \sQuote{modularity} object. }
}
\details{
  The input of this function should be a \sQuote{cna} class object
  containing \sQuote{network} and \sQuote{communities} attributes.

  This function reconstructs the community residue memberships for each
  modularity value. The purpose is to facilitate inspection of alternate 
  community partitioning points, which in practice often corresponds to a value 
  close to the maximum of the modularity, but not the maximum value itself.
}
\author{ Guido Scarabelli }  
\seealso{
   \code{\link{cna}}, \code{\link{network.amendment}}, \code{\link{summary.cna}}
}
\examples{
\donttest{
# PDB server connection required - testing excluded

if (!requireNamespace("igraph", quietly = TRUE)) {
   message('Need igraph installed to run this example')
} else {

###-- Build a CNA object
pdb <- read.pdb("4Q21")
modes <- nma(pdb)
cij <- dccm(modes)
net <- cna(cij, cutoff.cij=0.2)


##-- Reconstruct the community membership vector for each clustering step.
tree <- community.tree(net, rescale=TRUE)

## Plot modularity vs number of communities
plot( tree$num.of.comms, tree$modularity )

## Inspect the maximum modularity value partitioning
max.mod.ind <- which.max(tree$modularity)

## Number of communities (k) at max modularity
tree$num.of.comms[ max.mod.ind ]

## Membership vector at this partition point 
tree$tree[max.mod.ind,]

# Should be the same as that contained in the original CNA network object
net$communities$membership == tree$tree[max.mod.ind,]

# Inspect a new membership partitioning (at k=7)
memb.k7 <- tree$tree[ tree$num.of.comms == 7, ]

## Produce a new k=7 community network  
net.7 <- network.amendment(net, memb.k7)
plot(net.7, pdb)
#view.cna(net.7, trim.pdb(pdb, atom.select(pdb,"calpha")), launch=TRUE )

}
}
}
\keyword{analysis}
