\encoding{UTF-8}
\name{connectivity}
\alias{connectivity}
\alias{connectivity.default}
\alias{connectivity.pdb}
\title{Compute Connectivity From Coordinates}
\description{Search the atoms at a distance lower than the bonding distance (sum of covalent radii) to define the connectivity.}
\usage{
connectivity(...)

\method{connectivity}{default}(eleno.1, eleno.2, ...)

\method{connectivity}{connectivity}(x, ...)

\method{connectivity}{xyz}(x, ele.symb, safety = 1.2, by.block = FALSE, \dots)

\method{connectivity}{pdb}(x, atom.sel = atom.select(x, "notwater"), d.cut = 4, safety = 1.2, by.block = TRUE, \dots)
}
\arguments{
  \item{x}{an R object used to create a connectivity object.}
  \item{ele.symb}{a character vector containing atomic symbols used to find covalent radii.}
  \item{safety}{a numeric value used to extend the atomic radii.}
  \item{by.block}{a logical value indicating whether the connectivity has to be determine by block (see details).}
  \item{atom.sel}{an object of class \sQuote{select} indicating the atoms for which to compute the connectivity}
  \item{d.cut}{a numeric value. Use to compute Calpha connectivity (see \code{\link{calpha.connectivity}})}
  \item{\dots}{arguments passed to methods.}
}
\details{
The default method is a builder function which create an object of class \sQuote{connectivity} from its components (\code{eleno.1} and \code{eleno.2}). The method for object of class \sQuote{connectivity} is a validity function. It returns the same object as its argument if it is valid. Otherwise an error is returned.

For \sQuote{xyz} and \sQuote{pdb} objects, a distance matrix is computed, then, for each pair of atom the distance is compared to a bounding distance computed from atomic radii. If this distance is lower than the bounding distance then the atoms are assumed to be connected. If \code{by.block == TRUE}, a grid is defined to determined the connectivity by block. The method is slow but allow to deal with very large systems avoiding memory issues.

For \sQuote{pdb} objects, if only Caplha are found \code{calpha.connectivity} is called to compute the connectivity.
}
\value{
A two-columns \sQuote{connectivity} data.frame containing the indexes of bonded atoms. The column names are: \code{eleno.1} and \code{eleno.2}.
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
