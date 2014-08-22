\name{pca.3dalign}
\alias{pca.3dalign}
\title{ Principal Component Analysis }
\description{
   Performs principal components analysis (PCA) on an ensemble of PDB structures.
 }
\usage{
pca.3dalign(pdbs, core.find = FALSE, fit = FALSE, ...)
}
\arguments{
  \item{pdbs}{ an object of class \code{3dalign} as obtained from
    function \code{pdbaln} or \code{read.fasta.pdb}. }
  \item{core.find}{ logic, if TRUE core.find() function will be called to find core positions and coordinates of PDB structures will be fitted based on cores. }
  \item{fit}{ logic, if TRUE coordinates of PDB structures will be fitted based on all CA atoms. }
  \item{\dots}{ additional arguments passed to the method pca.xyz(). }
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
\seealso{ \code{\link{pca}}, \code{\link{pca.xyz}}}
\examples{

attach(transducin)

#-- Do PCA ignoring gap containing positions
pc.xray <- pca(pdbs, rm.gaps=TRUE)

# Plot results (conformer plots & scree plot)
plot(pc.xray, col=annotation[, "color"])

}
\keyword{ utilities }
\keyword{ multivariate }

