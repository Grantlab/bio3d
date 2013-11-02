\name{comms.significance}
\alias{comms.significance}
\title{
Wilcox Pair Test Calculation for Community Clustering
}
\description{
This function calculates the Wilcox Pair Test to determine community clustering goodness.
The input are the adjacency matrix, the community membership vector and (if any) the community membership vector names.
}
\usage{
comms.significance(adjmatrix, comms=comms, vnames=colnames(adjmatrix))
}
\arguments{
  \item{adjmatrix}{
Numeric adjacency matrix used to build the network
}
  \item{comms}{
Igraph community object. It must have a $membership attribute
}
  \item{vnames}{
(optional) Character vector containing the names of the network nodes.
}
}
\details{
This function calculates the Wilcox rank-sum test on the "internal" degree of a community (i.e. the edges within a community) 
and "external" degree of a community (i.e. the edges of a community with the rest of the network) 
to assess the quality of the community clustering.  The null hypothesis of the test is that there is no difference between
the two cases. The p-values of the test a printed in the output.
Be aware that a small sample size (i.e. small number of nodes in a community) might affecte the values calculated.
}
\value{
This function returns a numeric list with the p-value of Wilcox rank-sum test for each community.
}
\references{
More detailed instructions are available on this page:
http://igraph.wikidot.com/community-detection-in-r#toc2
}
\author{
Barry Grant, Guido Scarabelli
}
\examples{
## load input file
data(node.8)
pth <- system.file("data/", package = "cna")

## determine the community structure with igraph
igraph.comms.rslt <- igraph.comms(node.8)

## calculate the Wilcox test
wilcox.out <- comms.significance(node.8,igraph.comms.rslt)
}
\keyword{utility}

