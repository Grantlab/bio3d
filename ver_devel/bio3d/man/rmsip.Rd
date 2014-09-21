\name{rmsip}
\alias{rmsip}
\alias{rmsip.enma}
\alias{rmsip.nma}
\alias{rmsip.pca}
\alias{rmsip.pca.loadings}
\alias{rmsip.matrix}
\title{ Root Mean Square Inner Product }
\description{
  Calculate the RMSIP between two mode subspaces. 
}
\usage{
rmsip(...)
\method{rmsip}{nma}(...)
\method{rmsip}{pca}(...)
\method{rmsip}{pca.loadings}(...)
\method{rmsip}{enma}(enma, ncore=NULL, subset=10, ...)
\method{rmsip}{matrix}(modes.a, modes.b, subset=10,
                       row.name="a", col.name="b", ...)
}
\arguments{
  \item{enma}{ an object of class \code{"enma"} obtained from function
    \code{nma.pdbs}. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{subset}{ the number of modes to consider. }   
  \item{modes.a}{ an object of class \code{"pca"} or \code{"nma"} as
    obtained from functions \code{pca.xyz} or \code{nma}. }
  \item{modes.b}{ an object of class \code{"pca"} or \code{"nma"} as
    obtained from functions \code{pca.xyz} or \code{nma}. }
  \item{row.name}{ prefix name for the rows. }
  \item{col.name}{ prefix name for the columns. }
  \item{\dots}{ arguments passed to associated functions. }
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
  
  For function \code{rmsip.enma} a numeric matrix containing all
  pairwise RMSIP values of the modes stored in the \code{enma} object.
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
   Amadei, A. et al. (1999) \emph{Proteins} \bold{36}, 19--424.   
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{pca.xyz}}, \code{\link{nma}}, \code{\link{overlap}} }
\examples{
\dontrun{
# Load data for HIV example
trj <- read.dcd(system.file("examples/hivp.dcd", package="bio3d"))
pdb <- read.pdb(system.file("examples/hivp.pdb", package="bio3d"))

# Do PCA on simulation data
xyz.md <- fit.xyz(pdb$xyz, trj, fixed.inds=1:ncol(trj))
pc.sim <- pca.xyz(xyz.md)

# NMA 
modes <- nma(pdb)

# Calculate the RMSIP between the MD-PCs and the NMA-MODEs
r <- rmsip(modes, pc.sim, subset=10, row.name="NMA", col.name="PCA")

# Plot pairwise overlap values
plot(r, xlab="NMA", ylab="PCA")
}
}
\keyword{ utilities }
