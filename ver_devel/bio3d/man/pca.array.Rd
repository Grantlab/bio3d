\name{pca.array}
\alias{pca.array}
\title{ Principal Component Analysis of an arrary of matrices }
\description{
  Calculate the principal components of an array of correlation or
  covariance matrices.
}
\usage{
\method{pca}{array}(x)
}
\arguments{
  \item{x}{ an array of matrices, e.g. correlation or covariance
    matrices as obtained from functions \code{dccm} or \code{enma2covs}. }
}
\details{
  This function performs an PCA of residue-residue cross-correlations or
  covariance matrices derived from ensemble NMA of M structures. 
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
