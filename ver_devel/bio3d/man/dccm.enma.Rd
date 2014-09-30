\name{dccm.enma}
\alias{dccm.enma}
\title{ Cross-Correlation for Ensemble NMA (eNMA) }
\description{
  Calculate the cross-correlation matrices from an ensemble of NMA
  objects. 
}
\usage{
\method{dccm}{enma}(x, ncore = NULL, na.rm=FALSE, \dots)
}
\arguments{
  \item{x}{ an object of class \code{enma} as obtained from
    function \code{nma.pdbs}. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{na.rm }{ logical, if FALSE the DCCM might containt NA values
    (applies only when the \code{enma} object is calculated with
    argument \sQuote{rm.gaps=FALSE}). }
  \item{\dots}{ additional arguments passed to \code{dccm.nma}. }
}
\details{
  This is a wrapper function for calling \code{dccm.nma} on a collection
  of \sQuote{nma} objects as obtained from function \code{nma.pdbs}.
  
  See examples for more details.
}
\value{
  Returns a list with the following components:
  \item{all.dccm }{ an array or list containing the correlation matrices
    for each \sQuote{nma} object. An array is returned when the
    \sQuote{enma} object is calculated with \sQuote{rm.gaps=TRUE}, and a
    list is used when \sQuote{rm.gaps=FALSE}. }
  \item{avg.dccm }{ a numeric matrix containing the average correlation
    matrix. The average is only calculated when the \sQuote{enma} object
    is calculated with \sQuote{rm.gaps=TRUE}. }
}
\references{
  Wynsberghe. A.W.V, Cui, Q. \emph{Structure} \bold{14}, 1647--1653.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{ \code{\link{nma}}, \code{\link{dccm.nma}}, \code{\link{plot.dccm}} }
\examples{
## Fetch PDB files and split to chain A only PDB files
ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
raw.files <- get.pdb(ids, path = tempdir())
files <- pdbsplit(raw.files, ids, path = tempdir())

## Sequence/Structure Alignement
pdbs <- pdbaln(files)

## Normal mode analysis on aligned data
modes <- nma(pdbs)

## Calculate all 6 correlation matrices
cij <- dccm(modes)

## Plot correlations for first structure
plot.dccm(cij$all.dccm[,,1])

}
\keyword{ analysis }