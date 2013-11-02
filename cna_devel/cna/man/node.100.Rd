\name{node.100}
\alias{node.100}
\alias{node.50}
\alias{node.8}
\docType{data}
\title{
 Example Adjacency Matrices For The cna Package
}
\description{
A numeric nXn adjacency matrix, where n is the number of nodes. If two nodes are in contact, then the correspondig element in the matrix will be >0. If two nodes are not in contact the corresponding matrix element will be 0. The diagonal has been set to 0.
}
\usage{data(node.100)}
\format{
  The format is nXn numeric matrix
}
\examples{
data(node.100)
dim(node.100)

x<-cna(node.100)
plot(x)
plot(x$raw.network)
y=simplify(contract.vertices(x$raw.network, mapping=x$raw.communities$membership),
           edge.attr.comb=list(weight="max", name=toString, "ignore"))
plot(y)
}
\keyword{datasets}
