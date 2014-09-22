\name{pca}
\alias{pca}
\title{ Principal Component Analysis }
\description{
  Performs principal components analysis (PCA) on biomolecular structure data. 
}
\usage{
  pca(...)
}
\arguments{
  \item{\dots}{ arguments passed to the methods \code{pca.xyz}, \code{pca.pdbs}, etc. Typically this includes either a numeric matrix of Cartesian
    coordinates with a row per structure/frame (function
    \code{pca.xyz()}), or an object  of class \code{pdbs} as obtained from
    function \code{pdbaln} or \code{read.fasta.pdb} (function
    \code{pca.pdbs()}). } 
}
\details{
  Principal component analysis can be performed on any structure dataset of equal or unequal sequence composition to capture and characterize inter-conformer relationships.

  This generic \code{pca} function calls the corresponding methods function for actual calculation, which is determined by the class of the input argument \code{x}. Use
  \code{methods("pca")} to list all the current methods for \code{pca}
  generic. These will include:

  \code{\link{pca.xyz}}, which will be used when \code{x} is a numeric matrix
  containing Cartesian coordinates (e.g. trajectory data).
  
  \code{\link{pca.pdbs}}, which will perform PCA on the 
  Cartesian coordinates of a input \code{pdbs} object (as obtained from 
  the \sQuote{read.fasta.pdb} or \sQuote{pdbaln} functions).

  Currently, function \code{\link{pca.tor}} should be called explicitly as there
  are currently no defined \sQuote{tor} object classes.
  
  See the documentation and examples for each individual function for
  more details and worked examples.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant, Lars Skjaerven }
\seealso{
  \code{\link{pca.xyz}}, \code{\link{pca.pdbs}},
  \code{\link{pdbaln}}. }
\keyword{ utilities }
