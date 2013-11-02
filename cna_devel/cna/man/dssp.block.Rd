\name{dssp.block}
\alias{dssp.block}
\title{
  DSSP Block Generator to Collapse cij Matrices 
}
\description{
  This function creates a membership vector based on the secondary
  structure and collapse a cij matrix according to it.
}
\usage{
 dssp.block(cij=cij, sse=sse, collapse.method="max")
}
\arguments{
  \item{cij}{ A squared numeric matrix with NxN dimensions, where N is equal to
    the number of residues, containing cij couplings obtained with
    \sQuote{dccm} function. }
  \item{sse}{ Either a secondary structure object obtained with \sQuote{dssp}
    function or a numeric vector membership to be used in the collapsing
    process .}
  \item{collapse.method}{ A character object specifying the collapsing
    method to be used. The options available are: "max" (default),
    "median", "mean", "trimmed". }
}
\details{
  This function collapse a cij matrix according to the secondary
  structure of the protein.
}
\value{
 This function returns a list with the sse.membership vector and the
 cij.collapse matrix.
}
\references{
}
\author{ Guido and Barry }

\seealso{ \code{\link{dssp}}, \code{\link{dccm}} }
\examples{
library(cna)

## Load example files
pdb <- read.pdb("kin.pdb")
trj <- read.ncdf("kin.trj1.ncdf")

## Calculate secondary structure
sse <- dssp(pdb)

## Lets work with Calpha atoms only
ca.inds <- atom.select(pdb, "calpha")

## generation of the cij for each simulations (3 sims)
cij.1 <- dccm( trj[,ca.inds$xyz], ncore=8)

dssp.cij <- dssp.block(cij.1, sse.membership=sse)

plot.dccm(dssp.cij$cij.collapsed, at=c(-1,-0.8,-0.6,-0.4,0.4,0.6,0.8,1))

}
\keyword{ utility }