\name{rmsip}
\alias{rmsip}
\title{ Root Mean Square Inner Product }
\description{
  Calculate the RMSIP between two mode subspaces. 
}
\usage{
rmsip(modes.a, modes.b, subset=10, row.name="a", col.name="b")
}
\arguments{
  \item{modes.a}{ an object of class \code{"pca"} or \code{"nma"} as
    obtained from functions \code{pca.xyz} or \code{nma}. }
  \item{modes.b}{ an object of class \code{"pca"} or \code{"nma"} as
    obtained from functions \code{pca.xyz} or \code{nma}. }
  \item{subset}{ the number of modes to consider. }
  \item{row.name}{ prefix name for the rows. }
  \item{col.name}{ prefix name for the columns. }
}
\details{
  RMSIP is a measure for the similarity between two set of modes obtained
  from principal component or normal modes analysis. 
}
\value{
  Returns an \code{rmsip} object with the following components:
  \item{overlap}{ a numeric matrix containing pairwise (squared) dot products 
    between the modes. }
  \item{rmsip}{ a numeric RMSIP value. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
   Amadei, A. et al. (1999) \emph{Proteins} \bold{36}, 19--424.   
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{pca.xyz}}, \code{\link{nma}}, \code{\link{overlap}} }
\examples{
\dontrun{
# Do PCA on simulation data
pc.sim <- pca.xyz(xyz.md)

# NMA 
modes <- nma(pdb)

# Calculate the RMSIP between the Xray-PCs and the MD-PCs
r <- rmsip(modes, pc.sim, subset=10, row.name="NMA", col.name="PCA")

# Plot pairwise overlap values
plot(r, xlab="NMA", ylab="PCA")
}
}
\keyword{ utilities }
