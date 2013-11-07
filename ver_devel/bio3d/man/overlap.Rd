\name{overlap}
\alias{overlap}
\title{ Overlap analysis }
\description{
  Calculate the Squared Overlap between sets of vectors.
}
\usage{
overlap(modes, dv, nmodes=20)
}
\arguments{
  \item{modes}{ an object of class \code{"pca"} or \code{"nma"} as obtained from
    function \code{pca.xyz} or \code{nma}. Alternatively a 3NxM matrix
    of eigenvectors can be provided. } 
  \item{dv}{ a displacement vector of length 3N. } 
  \item{nmodes}{ the number of modes in which the calculation should be based. }
}
\details{
  Squared overlap (or dot product) is used to measure the similiarity 
  between a displacement vector (e.g. a difference vector between two 
  conformational states) and mode vectors obtained from principal component 
  or normal modes analysis. 

  By definition the cumulative sum of the overlap values equals to one. 

  Structure \code{modes$U} (or alternatively, the 3NxM matrix of eigenvectors)
  should be of same length (3N) as \code{dv}. 
}
\value{
  Returns a list with the following components:
  \item{overlap}{ a numeric vector of the squared dot products (overlap values)
    between the (normalized) vector (\code{dv}) and each mode in \code{mode}. }
  \item{overlap.cum}{ a numeric vector of the cumulative squared overlap values. }
}
\references{
  Skjaerven, L. et al. (2011) \emph{Proteins} \bold{79}, 232--243.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{rmsip}}, \code{\link{pca.xyz}}, \code{\link{nma}},
  \code{\link{difference.vector}} }
\examples{
data(kinesin)
attach(kinesin, warn.conflicts=FALSE)
     
# Ignore gap containing positions
##gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)
     
#-- Do PCA
pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])
    
# Define a difference vector between two structural states
diff.inds <- c(grep("d1v8ja", pdbs$id), 
               grep("d1goja", pdbs$id))

dv <- difference.vector( pdbs$xyz[diff.inds,], gaps.pos$f.inds )

# Calculate the squared overlap between the PCs and the difference vector  
o <- overlap(pc.xray, dv)
o <- overlap(pc.xray$U, dv)


## Calculate overlap from NMA
pdb.a <- read.pdb("1cmk")
pdb.b <- read.pdb("3dnd")

## Fetch CA coordinates
sele.a <- atom.select(pdb.a, "//E/15:350///CA")
sele.b <- atom.select(pdb.b, "//A/1:350///CA")

xyz <- rbind(pdb.a$xyz[sele.a$xyz],
             pdb.b$xyz[sele.b$xyz])

## Superimpose
xyz[2,] <- fit.xyz(xyz[1,], xyz[2,], 1:ncol(xyz))

## The difference between the two conformations
dv <- difference.vector( xyz )

## Calculate normal modes
modes <- nma(pdb.a, inds=sele.a)

# Calculate the squared overlap between the normal modes 
# and the difference vector  
o <- overlap(modes, dv)

detach(kinesin)
}
\keyword{ utilities }
