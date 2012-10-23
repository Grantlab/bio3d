\name{overlap}
\alias{overlap}
\title{ Overlap analysis }
\description{
  Calculate the Squared Overlap between sets of vectors
}
\usage{
overlap(pca, dv, num.modes=20)
}
\arguments{
  \item{pca}{ an object of class \code{"pca"} as obtained from
    function \code{pca.xyz}, or alternatively a 3NxM matrix of eigenvectors. } 
  \item{dv}{ a displacement vector of length 3N. } 
  \item{num.modes}{ the number of modes to consider. }
}
\details{
  Squared overlap (or dot product) is used to measure the similiarity 
  between a displacement vector (e.g. a difference vector between two 
  conformational states) and mode vectors obtained from principal component 
  or normal modes analysis. 

  By definition the cumulative sum of the overlap values equals to one. 

  Structure \code{pca$U} (or alternatively, the 3NxM matrix of eigenvectors)
  should be of same length (3N) as \code{dv}. 
}
\value{
  Returns a list with the following components:
  \item{overlap}{ a numeric vector of the squared dot products (overlap values)
    between the (normalized) vector (\code{dv}) and each mode in \code{pca}. }
  \item{overlap.cum}{ a numeric vector of the cumulative squared overlap values. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
   }
\author{ Lars Skjærven }
\seealso{ \code{\link{rmsip}}, \code{\link{pca.xyz}} }
\examples{
     data(kinesin)
     attach(kinesin)
     
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
     o <- overlap(pc.xray$U, dv)
}
\keyword{ utilities }
