\name{rmsip}
\alias{rmsip}
\title{ Root Mean Square Inner Product }
\description{
  Calculate the RMSIP between two mode subspaces. 
}
\usage{
rmsip(pca.a, pca.b, subset=10, row.name=NULL, col.name=NULL)
}
\arguments{
  \item{pca.a}{ an object of class \code{"pca"} as obtained from
    functions \code{pca.xyz}. }
  \item{pca.b}{ an object of class \code{"pca"} as obtained from
    functions \code{pca.xyz}. }
  \item{subset}{ the number of modes to consider. }
  \item{row.name}{ prefix name for the rows. }
  \item{col.name}{ prefix name for the columns. }
}
\details{
  RMSIP is a measure for the similarity between two set of modes obtained
  from principal component or normal modes analysis. 

  Structure \code{pca.a$U} and \code{pca.b$U} should have the same number of columns.
   
}
\value{
  Returns a list with the following components:
  \item{overlap}{ a numeric matrix containing pairwise (squared) dot products 
    between the modes. }
  \item{rmsip}{ a numeric RMSIP value. }
}
\references{
   Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
   Amadei, A. et al. (1999) \emph{Proteins} \bold{36}, 19--424.   
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{overlap}}, \code{\link{pca.xyz}} }
\examples{
     data(kinesin)
     attach(kinesin, warn.conflicts=FALSE)
     
     # Ignore gap containing positions
     gaps.pos <- gap.inspect(pdbs$xyz)
     
     # Do PCA
     pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
     
     \dontrun{
     # Read in MD trajectory
     sim <- NULL
     sim$trj <- read.ncdf("300ns_CA.nc")
     sim$xyz <- fit.xyz(fixed = pdbs$xyz[1, ], mobile =  sim$trj,
     	                fixed.inds = core$c1A.xyz, mobile.inds = core$c1A.xyz,
                        full.pdbs = FALSE, het2atom = TRUE)

     # Do PCA on the simulation
     sim$pc <- pca.xyz(sim$xyz)

     # Overlap analysis
     o <- overlap(sim$pc, dv, num.modes=50)

     # Calculate the RMSIP between the Xray-PCs and the MD-PCs
     r <- rmsip(pc.xray, sim$pc, subset=10, row.name="Xray", col.name="MD")

     # Plot pairwise overlap values 
     image(1:10, 1:10, r$overlap, col=gray(50:0/50), zlim=c(0,1),
           xlab="MD-PCs", ylab="X-ray PCs")

     # Convergence of eigenvectors
     all.rmsip <- c()
     intervals <- seq(100,1000, by=100)
     for ( i in 1:length(intervals) ) {
        traj.inds <- 1:intervals[i]
        tmp.pc <- pca.xyz(sim$xyz[traj.inds,])
        r <- rmsip(sim$pc, tmp.pc, subset=10)
        all.rmsip <- c(all.rmsip, r$rmsip)
     }

}
}
\keyword{ utilities }
