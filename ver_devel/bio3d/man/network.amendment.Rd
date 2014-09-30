\name{network.amendment}
\alias{network.amendment}
\title{ Amendment of a CNA Network According To A Input Community Membership Vector.  }
\description{
  This function changes the \sQuote{communities} attribute of a 
  \sQuote{cna} class object to match a given membership vector.
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
  This function is useful, in combination with \sQuote{community.tree}, 
  for inspecting different community partitioning options of a input 
  \sQuote{cna} object. See examples.
} 
\author{ Guido Scarabelli }  
\seealso{
   \code{\link{cna}}, \code{\link{community.tree}}, \code{\link{summary.cna}} 
}
\examples{

##-- Build a CNA object
pdb <- read.pdb("4Q21")
modes <- nma(pdb)
cij <- dccm(modes)
net <- cna(cij, cutoff.cij=0.2)

##-- Community membership vector for each clustering step
tree <- community.tree(net, rescale=TRUE)

## Produce a new k=7 membership vector and CNA network
memb.k7 <- tree$tree[ tree$num.of.comms == 7, ]
net.7 <- network.amendment(net, memb.k7)

plot(net.7, pdb)

print(net)
print(net.7)

}
\keyword{utility}
