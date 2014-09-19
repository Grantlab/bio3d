\name{filter.dccm}
\alias{filter.dccm}
\title{
  Filter for Cross-correlation Matrices (Cij)
}
\description{
  This function builds various cij matrix for dynamical network analysis
}
\usage{
  filtr.dccm(x, xyz, fac = NULL, cutoff.cij = 0.4, collapse = TRUE,  
      extra.filter = NULL, ...)
}
\arguments{
  \item{x}{ A matrix (nXn), a numeric array with 3 dimensions (nXnXm), 
    a list with m cells each containing nXn matrix, or a list with
    \sQuote{all.dccm} component, containing atomic correlation values, 
    where "n" is the residue number and "m" the system/frame number. 
    The matrix elements should be in between -1 and 1. See \sQuote{dccm} 
    function in bio3d package for further details. }
  \item{xyz}{ XYZ coordinates for distance matrix calculation. }
  \item{fac}{ factor indicating distinct categories of input correlation matrices. }
  \item{cutoff.cij}{ Threshold to define significant couplings. }
  \item{collapse}{ logical, if TRUE the mean matrix will be returned. }
  \item{extra.filter}{ Filter to apply in addition to the model chosen. } 
  \item{\dots}{ extra arguments passed to function \code{cmap}. } 
}
\value{
  Returns a matrix or a 3D array of filtered cij.
}
\details{ 
     The function inspects a set of correlation matrices and decides edges for 
     correlation network analysis based on:
     1. min(abs(cij)) >= cutoff.cij, or
     2. max(abs(cij)) >= cutoff.cij && residues contact each other based on results
     from \code{cmap}.
}  
\seealso{
   \code{\link{cna}}, \code{\link{dccm}}, \code{\link{cmap}} 
}
\examples{

\dontrun{
data(transducin)
attach(transducin)

gaps.pos <- gap.inspect(pdbs$xyz)
modes <- nma.pdbs(pdbs, full=TRUE)
dccms <- dccm.enma(modes)

cij <- filter.dccm(dccms, xyz=pdbs)

}
}
\keyword{analysis}
