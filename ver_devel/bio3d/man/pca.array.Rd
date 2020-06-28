\name{pca.array}
\alias{pca.array}
\title{ Principal Component Analysis of an array of matrices }
\description{
  Calculate the principal components of an array of correlation or
  covariance matrices.
}
\usage{
\method{pca}{array}(x, use.svd = TRUE, rm.gaps=TRUE, ...)
}
\arguments{
  \item{x}{ an array of matrices, e.g. correlation or covariance
    matrices as obtained from functions \code{dccm} or \code{enma2covs}. }
  \item{use.svd}{ logical, if TRUE singular value decomposition (SVD) is 
    called instead of eigenvalue decomposition. }
  \item{rm.gaps}{ logical, if TRUE gap cells (with missing
    coordinate data in any input matrix) are removed before
    calculation. This is equivalent to removing NA cells from x. }
  \item{\dots}{ . }
}
\details{
  This function performs PCA of symmetric matrices, such as distance 
  matrices from an ensemble of crystallographic structures, residue-residue 
  cross-correlations or covariance matrices derived from ensemble NMA 
  or MD simulation replicates, and so on. The \sQuote{upper triangular} 
  region of the matrix is regarded as a long vector of random variables. 
  The function returns M eigenvalues and eigenvectors with each eigenvector 
  having the dimension N(N-1)/2, where M is the number of matrices and N 
  the number of rows/columns of matrices.
}
\value{
  Returns a list with components equivalent to the output from
  \code{pca.xyz}. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.   
}
\author{ Xin-Qiu Yao, Lars Skjaerven }
\seealso{ \code{\link{pca.xyz}} }
\keyword{ utilities }
