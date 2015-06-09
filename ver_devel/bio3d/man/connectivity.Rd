\encoding{UTF-8}
\name{connectivity}
\alias{are.symb}
\alias{connectivity}
\alias{connectivity.default}
\alias{connectivity.pdb}
\title{Compute Connectivity From Coordinates}
\description{Search the atoms at a distance lower than the bonding distance (sum of covalent radii) to define the connectivity.}
\usage{
connectivity(...)

\method{connectivity}{default}(eleno.1, eleno.2, ...)

\method{connectivity}{connectivity}(x, ...)

\method{connectivity}{xyz}(x, elesy, safety = 1.2, by.block = FALSE, ca.check = FALSE, \dots)

\method{connectivity}{pdb}(x, safety = 1.2, by.block = TRUE, \dots)

connectivity(object) <- value

\method{connectivity}{pdb}(object) <- value

are.symb(x)

}
\arguments{
  \item{x}{an R object, typically a \sQuote{pdb} object or \sQuote{xyz} coordinate vector for which connectivity is to be determined or set.}
  \item{elesy}{a character vector containing atomic symbols used to find covalent radii.}
  \item{safety}{a numeric value used to extend the atomic radii.}
  \item{by.block}{a logical value indicating whether the connectivity has to be determine by block (see details).}
  \item{ca.check}{a logical value indicating whether to check for C-alpha atom connectivity based on supplied elesy. Used primarily for visualization purposes (see \code{\link{calpha.connectivity}})} 
  \item{object}{an R object, typically a \sQuote{pdb} object, whose \sQuote{connectivity} component will be replaced.}
  \item{value}{an object of class \sQuote{connectivity}.}
  \item{eleno.1}{a vector of names for connectivity elements. }
  \item{eleno.2}{a vector of names for connectivity elements. }
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
