\name{pca.pdbs}
\alias{pca.pdbs}
\title{ Principal Component Analysis }
\description{
  Performs principal components analysis (PCA) on an ensemble of PDB structures.
}
\usage{
\method{pca}{pdbs}(pdbs, core.find = FALSE, fit = FALSE, \dots)
}
\arguments{
  \item{pdbs}{ an object of class \code{pdbs} as obtained from
    function \code{pdbaln} or \code{read.fasta.pdb}. }
  \item{core.find}{ logical, if TRUE core.find() function will be called
    to find core positions and coordinates of PDB structures will be
    fitted based on cores.  }
  \item{fit}{ logical, if TRUE coordinates of PDB structures will be
    fitted based on all CA atoms.  }
  \item{\dots}{ additional arguments passed to the method \code{pca.xyz}. }
}
\details{
  The function \code{pca.pdbs} is a wrapper for the function
  \code{\link{pca.xyz}}, wherein more details of the PCA procedure
  are documented. 
}
\value{
  Returns a list with the following components:
  \item{L }{eigenvalues.}
  \item{U }{eigenvectors (i.e. the variable loadings).}
  \item{z.u }{scores of the supplied \code{data} on the pcs.}
  \item{sdev }{the standard deviations of the pcs.}
  \item{mean }{the means that were subtracted.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{  Barry Grant and Karim ElSawy }
\seealso{
  \code{\link{pca}}, \code{\link{pca.xyz}},
  \code{\link{pdbaln}}, \code{\link{nma}}. 
}
\examples{

attach(transducin)

#-- Do PCA ignoring gap containing positions
pc.xray <- pca(pdbs)

# Plot results (conformer plots & scree plot)
plot(pc.xray, col=annotation[, "color"])

}
\keyword{ utilities }
\keyword{ multivariate }

