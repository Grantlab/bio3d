\name{cij.filter}
\alias{cij.filter}
\title{
  Filter for Cross-correlation Matrices (Cij)
}
\description{
  This function builds various cij matrix for dynamical network analysis
}
\usage{
  cij.filter(cij, inds = 1:dim(cij)[3L], xyz = NULL, 
      model = c("simple", "minimal", "full", "pca"), 
      nmodes = 3, cutoff.var = 0.0025, cutoff.cij = 0.3, 
      cutoff.loading = 0.02, cutoff.dm = 15, extra.filter = NULL)
}
\arguments{
  \item{cij}{ A matrix (nXn), a numeric array with 3 dimensions (nXnXm), 
    a list with m cells each containing nXn matrix, or a list with
    \sQuote{all.dccm} component, containing atomic correlation values, 
    where "n" is the residue number and "m" the system/frame number. 
    The matrix elements should be in between -1 and 1. See \sQuote{dccm} 
    function in bio3d package for further details. }
  \item{inds}{ Indices for a subset to operate on. }
  \item{xyz}{ XYZ coordinates for distance matrix calculation. }
  \item{model}{ Model to choose for filtering. }
  \item{nmodes}{ Number of modes to add. Used for model = "pca". }
  \item{cutoff.var}{ Threshold to define variable couplings. }
  \item{cutoff.cij}{ Threshold to define significant couplings. }
  \item{cutoff.loading}{ Threshold to define significant couplings. 
                  Used for model = "pca". }
  \item{cutoff.dm}{ Threshold to define structure contact. } 
  \item{extra.filter}{ Filter to apply in addition to the model chosen. } 
}
\value{
  Returns a matrix or a 3D array of filtered cij.
}
\details{ 
  MODELS:
  \item{simple}{ Remove edges with abs(cij) < cutoff.cij. }
  \item{minimal}{ 1. Remove edges with max(abs(cij)) < cutoff.cij; 
     2. Remove edges with min(abs(cij)) < cutoff.cij && min(dij) > cutoff.dm; 
     3. Remove edges with var(abs(cij)) >= cutoff.var. }
  \item{full}{ 1. Remove edges with max(abs(cij)) < cutoff.cij;
     2. Remove edges with min(abs(cij)) < cutoff.cij && min(dij) > cutoff.dm. }
  \item{pca}{ 1. Remove edges with max(abs(cij)) < cutoff.cij; 
     2. Remove edges with min(abs(cij)) < cutoff.cij && min(dij) > cutoff.dm; 
     3. Remove edges with var(abs(cij)) >= cutoff.var && intersect({pc-loading_n(ij) < cutoff.loading; n=1:nmodes}). }
}  
\seealso{
   \code{\link{cna}}, \code{\link{dccm}}
}
\examples{

\dontrun{
data(transducin)
attach(transducin)

gaps.pos <- gap.inspect(pdbs$xyz)
modes <- nma.pdbs(pdbs, full=TRUE)
dccms <- dccm.enma(modes)

cij1 <- cij.filter(dccms, inds = 1, model = "simple")
cij2 <- cij.filter(dccms, xyz = pdbs, model = "minimal")
cij3 <- cij.filter(dccms, xyz = pdbs, model = "full")
cij4 <- cij.filter(dccms, xyz = pdbs, model = "pca", nmodes = 3)

}
}
\keyword{analysis}
