\name{difference.vector}
\alias{difference.vector}
\title{ Difference Vector }
\description{
  Define a difference vector between two conformational states.
}
\usage{
difference.vector(xyz, xyz.inds=NULL)
}
\arguments{
  \item{xyz}{ numeric matrix of Cartesian coordinates with a row per
      structure. } 
  \item{xyz.inds}{ a vector of indices that selects the elements of columns upon
          which the calculation should be based. }
}
\details{
  Squared overlap (or dot product) is used to measure the similiarity 
  between a displacement vector (e.g. a difference vector between two 
  conformational states) and mode vectors obtained from principal component 
  or normal modes analysis. 
}
\value{
  Returns a numeric vector (normalized) of the structural difference.
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
   }
\author{ Lars Skjaerven }
\seealso{ \code{\link{overlap}} }
\examples{
     data(kinesin)
     attach(kinesin, warn.conflicts=FALSE)
     
     # Ignore gap containing positions
     gaps.res <- gap.inspect(pdbs$ali)
     gaps.pos <- gap.inspect(pdbs$xyz)
     
     #-- Do PCA
     pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
     
     # Plot results (conformer plots & scree plot)
     plot(pc.xray)

    # Define a difference vector between two structural states
    diff.inds <- c(grep("d1v8ja", pdbs$id), 
                   grep("d1goja", pdbs$id))

    dv <- difference.vector( xyz[diff.inds,], gaps.pos$f.inds )

    # Calculate the squared overlap between the PCs and the difference vector  
    o <- overlap(pc.xray, dv)
}
\keyword{ utilities }
