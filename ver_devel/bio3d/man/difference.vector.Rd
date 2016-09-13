\name{difference.vector}
\alias{difference.vector}
\title{ Difference Vector }
\description{
  Define a difference vector between two conformational states.
}
\usage{
difference.vector(xyz, xyz.inds=NULL, normalize=FALSE)
}
\arguments{
  \item{xyz}{ numeric matrix of Cartesian coordinates with a row per
      structure. } 
  \item{xyz.inds}{ a vector of indices that selects the elements of columns upon
    which the calculation should be based. }
  \item{normalize}{ logical, if TRUE the difference vector is
    normalized. }
}
\details{
  Squared overlap (or dot product) is used to measure the similiarity 
  between a displacement vector (e.g. a difference vector between two 
  conformational states) and mode vectors obtained from principal component 
  or normal modes analysis. 
}
\value{
  Returns a numeric vector of the structural difference (normalized if desired).
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{overlap}} }
\examples{
attach(kinesin)
     
# Ignore gap containing positions
gaps.pos <- gap.inspect(pdbs$xyz)
     
#-- Do PCA
pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])
     
# Define a difference vector between two structural states
diff.inds <- c(grep("d1v8ja", pdbs$id), 
               grep("d1goja", pdbs$id))

## Calculate the difference vector
dv <- difference.vector( pdbs$xyz[diff.inds,], gaps.pos$f.inds )

# Calculate the squared overlap between the PCs and the difference vector  
o <- overlap(pc.xray, dv)

detach(kinesin)

}
\keyword{ utilities }
