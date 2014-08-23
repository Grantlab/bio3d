\name{nma}
\alias{nma}
\title{ Normal Mode Analysis }
\description{
  Performs normal mode analysis (NMA) on either a single or an ensemble
  of protein structures
}
\usage{
nma(x, ...)
}
\arguments{
  \item{x}{ an object of class \code{pdb} as obtained from
    function \code{read.pdb}. Alternatively, an object of class
    \code{pdbs} as obtained from function \code{pdbaln} or
    \code{read.fasta.pdb} that will be passed to the
    \code{nma.pdbs()} function, see below for examples. }
  \item{\dots}{ additional arguments passed to the methods
    \code{nma.pdb}, or \code{nma.pdbs}. }
}
\details{
  \code{nma} is a generic function calling the corresponding function
  determined by the class of the input argument \code{x}. Use
  \code{methods("nma")} to get all the methods for \code{nma}
  generic:

  \code{\link{nma.pdb}} will be used when \code{x} is of class
  \code{pdb} as obtained from function \code{read.pdb}.

  \code{\link{nma.pdbs}} will perform \sQuote{ensemble} normal mode
  analysis of the PDB structures stored in the \code{pdbs} object.
  
  See documentation and examples for each corresponding function for
  more details. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{nma.pdb}}, \code{\link{nma.pdbs}},
  \code{\link{pca}}.
}
\keyword{ analysis }
