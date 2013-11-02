\name{vmd.btwn}
\alias{vmd.btwn}
\title{
VMD edge betweenness vector
}
\description{
This function reads a VMD-DNA-plugin betweenness output files and creates a numeric edge betweenness vector.
}
\usage{
vmd.btwn(btwn.matrix = btwn.matrix)
}
\arguments{
  \item{btwn.matrix}{
Edge betweenness matrix output file from VMD-DNA-plugin. The matrix is square nXn, where "n" is the number of residues.
The values > 0 represent the edge betweenness for the residue connections. Values equal to 0 correspond to no connection.
Note: the first line of the file (head) is discarded.
}
}
\value{
This function returns a numeric vector with the edge betweenness values.
}

\examples{
## load input
pth <- system.file("data",package="cna")  ## path to the example files
vmd.btwn.8.node.input <- data.matrix(read.table(paste(pth,"/vmd-btwn-8-node.dat",sep=""),skip=1))  ## read vmd btwn matrix

## edge betweenness vector creation
vmd.btwn.rslt <- vmd.btwn(vmd.btwn.8.node.input) 

}
\keyword{utility}

