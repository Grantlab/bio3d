\name{igraph.btwn}
\alias{igraph.btwn}
\title{
Igraph Package Edge Betweenness
}
\description{
This function calculates the edge betweenness values using the Brandes algorithm (igraph package).
The higher is the betweenness value of an edge the more crucial is that edge for inter-community connection.
}
\usage{
igraph.btwn(network.matrix = network.matrix)
}
\arguments{
  \item{network.matrix}{
a numeric adjacency matrix object. i.e. the matrix elements are >0 if two elements are in contact, 0 if they are not.
}
}
\value{
This function returns a numeric edge betweenness vector. Each element is the betweenness value of an edge.
To show the edge list, use E(network) command.
}

\seealso{
graph.adjacency, edge.betweenness
}
\examples{
## load input
data(node.8) 

## calculate the edge betweenness vector
igraph.btwn.rslt <- igraph.btwn(node.8)

}
\keyword{analysis}

