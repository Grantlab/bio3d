\name{calpha.connectivity}
\alias{calpha.connectivity}
\alias{calpha.connectivity.xyz}
\alias{calpha.connectivity.pdb}
\title{Compute Calpha Atom Connectivity From Coordinates}
\description{Connects consecutive Calpha atoms at a distance lower than a given cutoff distance}
\usage{
calpha.connectivity(...)
\method{calpha.connectivity}{xyz}(x, d.cut = 4, ...)
\method{calpha.connectivity}{pdb}(x, d.cut = 4, ...)
}
\arguments{
  \item{x}{an R object for which to compute the connectivity.}
  \item{d.cut}{a numeric value. The cutoff distance after which atoms are not connected}
  \item{\dots}{arguments passed to methods.}
}
\details{
This function is primarily used for visualization purposes in the \sQuote{view} family of functions when dealing with C-alpha only structures or rendering chain traces. 

Operation wise, \code{calpha.connectivity} first calls \code{dist.xyz} to compute inter-atomic distances between consecutive atoms and then determines if they are connected (d < d.cut). Increasing d.cut will retain more connections when viewing interpolations and displacements from PCA and NMA for example.
}
\value{
A two-columns \sQuote{connectivity} data.frame containing the indexes of bonded atoms. The column names are: \code{eleno.1} and \code{eleno.2}.
}

\author{Julien Ide and Barry Grant}

\seealso{
\code{\link{connectivity}}
}
\examples{
x <- read.pdb("1bg2")
con <- calpha.connectivity(x)
}
\keyword{ manip }
