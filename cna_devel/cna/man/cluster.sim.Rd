\name{cluster.sim}
\alias{cluster.sim}
\title{
Cluster Similarity
}
\description{
This function calculates Jaccard, Rand and Adjusted Rand index for clustering similarity.
The input are two named vectors corresponding to the clustering membership.
}
\usage{
cluster.sim(grp1, grp2, self = FALSE, minSZ = 1, method = "jaccard")
}
\arguments{
  \item{grp1}{
named vector where the cluster IDs are values and the item labels are names
}
  \item{grp2}{
named vector where the cluster IDs are values and the item labels are names
}
  \item{self}{
FALSE to ignore or TRUE to include clusters with single items (default=FALSE)
}
  \item{minSZ}{
selection of a minimum cluster size (default=1)
}
  \item{method}{
supported similarity methods are "jaccard", "rand" and "arand" (default="jaccard")
}
}
\details{
Definition of Indices:\cr
- Jaccard Index: J(C,K) = a/(a+b+c)\cr
- Rand Index: R(C,K) = (a+d)/(a+b+c+d)\cr
- Adjusted Rand Index: AR(C,K) = (2*(a*d-c*b)) / ((a+b)*(b+d)+(a+c)*(c+d))\cr 
Where:\cr
a = number of pairs that cluster together in both C and K\cr
b = number of pairs that cluster together in C, but not K\cr
c = number of pairs that cluster together in K, but not C\cr
d = number of pairs that are not joined into clusters in K and C
}
\value{
This function returns a numeric list with:\cr
- The overlapping matrix ($intersects)\cr
- The variables used for the index calculation ($variables)\cr
- The index calculated (either $jaccard, $rand or $arand)\cr
}
\references{
More detailed instructions are available on this page:
http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#clustering_jaccard
}
\author{
Barry Grant, Guido Scarabelli
}
\examples{
\dontrun{
## load input file
data(node.8)
pth <- system.file("data/", package = "cna")

## determine the community structure with igraph
igraph.comms.rslt <- igraph.comms(node.8)

## read vmd community list
vmd.comms.rslt <- vmd.comms(paste(pth,"/vmd-8-node.out",sep=""),node.num=8)

## calcuate Jaccard index
jcc.index <- cluster.sim(igraph.comms.rslt$membership, vmd.comms.rslt$membership, self = FALSE, minSZ = 1, method = "jaccard")
}
}
\keyword{utility}

