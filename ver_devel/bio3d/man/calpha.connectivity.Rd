\name{calpha.connectivity}
\alias{calpha.connectivity}
\alias{calpha.connectivity.xyz}
\alias{calpha.connectivit.pdb}
\title{Compute Calpha Connectivity From Coordinates}
\description{Connects consecutive Calpha atoms at a distance lower than a given cutoff distance}
\usage{
calpha.connectivity(...)
\method{calpha.connectivity}{xyz}(x, d.cut = 4, ...)
\method{calpha.connectivity}{pdb}(x, atom.sel = atom.select(x, "protein"), d.cut = 4, ...)
}
\arguments{
  \item{x}{an R object for which to compute the connectivity.}
    \item{atom.sel}{an object of class \sQuote{select} indicating the atoms for which to compute the connectivity}
  \item{d.cut}{a numeric value. The cutoff distance after which atoms are not connected}
  \item{\dots}{arguments passed to methods.}
}
\details{
\code{calpha.connectivity} first call \code{trim.pdb} to extract Calpha. Then, \code{dist.xyz} is called to compute inter-atomic distances between consecutive atoms and determine if they are connected (d < d.cut).
}
\value{
A two-columns \sQuote{connectivity} data.frame containing the indexes of bonded atoms. The column names are: \code{eleno.1} and \code{eleno.2}.
}

\author{Julien Idé}

\seealso{
\code{\link{connectivity}}
}
\examples{
x <- read.pdb("1bg2")
con <- calpha.connectivity(x)
}
\keyword{ manip }
