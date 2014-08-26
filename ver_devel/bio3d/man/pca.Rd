\name{pca}
\alias{pca}
\title{ Principal Component Analysis }
\description{
  Performs principal components analysis (PCA) on protein structure data. 
}
\usage{
pca(...)
}
\arguments{
  \item{\dots}{ arguments passed to the methods \code{pca.xyz} or
    \code{pca.pdbs}. Provide either a numeric matrix of Cartesian
    coordinates with a row per structure/frame (function
    \code{pca.xyz()}), or an object  of class \code{pdbs} as obtained from
    function \code{pdbaln} or \code{read.fasta.pdb} (function
    \code{pca.pdbs()}). } 
}
\details{
  \code{pca} is a generic function calling the corresponding function
  determined by the class of the input argument \code{x}. Use
  \code{methods("pca")} to get all the methods for \code{pca}
  generic:

  \code{\link{pca.xyz}} will be used when \code{x} is a numeric matrix
  containing Cartesian coordinates (e.g. trajectory data).
  
  \code{\link{pca.pdbs}} will calculate the principal component
  of the Cartesian coordinates of the an \code{pdbs} object.

  Function \code{\link{pca.tor}} should be called explicitly since there
  are currently no \sQuote{tor} object.
  
  See documentation and examples for each corresponding function for
  more details.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{
  \code{\link{pca.xyz}}, \code{\link{pca.pdbs}},
  \code{\link{pdbaln}}. }
\keyword{ utilities }
