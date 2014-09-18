\encoding{UTF-8}
\name{connectivity}
\alias{connectivity}
\alias{connectivity.default}
\alias{connectivity.pdb}
\title{Compute Connectivity From Coordinates}
\description{Search the atoms at a distance lower than the bonding distance (sum of covalent radii) to define the connectivity.}
\usage{
connectivity(...)

\method{connectivity}{default}(x, ele.symb, safety = 1.2, by.block = FALSE, \dots)

\method{connectivity}{pdb}(x, safety = 1.2, by.block = FALSE, \dots)
}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{ele.symb}{a character vector containing atomic symbols used to find covalent radii.}
  \item{safety}{a numeric value used to extend the atomic radii.}
  \item{by.block}{a logical value indicating whether the connectivity has to be determine by block (see details).}
  \item{\dots}{arguments passed to methods.}
}
\details{
A distance matrix is computed, then, for each pair of atom the distance is compared to a bounding distance computed from atomic radii. If this distance is lower than the bounding distance then the atoms are assumed to be connected. If \code{by.block == TRUE}, a grid is defined to determined the connectivity by block. The method is slow but allow to deal with very large systems avoiding memory issues. 
}
\value{
A two-columns data.frame containing the indexes of bonded atoms. The column names are: \code{eleno.1} and \code{eleno.2}.
}

\author{Julien Idé}

\seealso{
\code{\link{elements}}, \code{\link{pdb}}
}
\examples{
x <- read.pdb("1bg2")
system.time(c1 <- connectivity(x))
system.time(c2 <- connectivity(x, by.block = TRUE))
all(c1==c2)
}
\keyword{ manip }
