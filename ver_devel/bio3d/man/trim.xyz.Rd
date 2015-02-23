\name{trim.xyz}
\alias{trim.xyz}
\title{ Trim a XYZ Object of Cartesian Coordinates. }
\description{
  Produce a new smaller XYZ object, containing a subset of atoms.
}
\usage{
\method{trim}{xyz}(xyz, row.inds = NULL, col.inds = NULL, \dots)
}
\arguments{
  \item{xyz}{ a XYZ object containing Cartesian coordinates,
    e.g. obtained from \code{\link{read.pdb}}, \code{\link{read.ncdf}}. }    
  \item{row.inds}{ a numeric vector specifying which rows of the xyz
    matrix to return. }
  \item{col.inds}{ a numeric vector specifying which columns of the xyz
    matrix to return. }
  \item{\dots}{ additional arguments passed to and from functions. } 
}
\details{
  This function provides basic functionality for subsetting a
  matrix of class \sQuote{xyz} while also maintaining the class
  attribute. 
}
\value{
  Returns an object of class \code{xyz} with the Cartesian coordinates
  stored in a matrix object with dimensions M x 3N, where N is the
  number of atoms, and M number of frames. 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{read.pdb}}, \code{\link{as.xyz}}.
}
\examples{
\dontrun{
## Read a PDB file from the RCSB online database
pdb <- read.pdb("1bg2")

## Select calpha atoms
sele <- atom.select(pdb, "calpha")

## Trim XYZ
trim(pdb$xyz, col.inds=sele$xyz)

## Equals to
pdb$xyz[, sele$xyz, drop=FALSE]
}
}
\keyword{ utilities }
