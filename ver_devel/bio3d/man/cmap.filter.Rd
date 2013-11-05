\name{cmap.filter}
\alias{cmap.filter}
\title{
Contact Map Consensus Filtering
}
\description{
This function filters a tridimensional contact matrix (nXnXz, where n is the residue number and z is the simulation number) selecting only the contact present in at least p simulations, where p<=z.
}
\usage{
cmap.filter(cm, cutoff.sims = dim(cm)[3])
}
\arguments{
  \item{cm}{ A numeric array with 3 dimensions (nXnXz) containing binary contact values. "n" is the residue number, "z" the simulation number. The matrix elements should be 1 if two residues are in contact and 0 if they are not in contact. 
}
  \item{cutoff.sims}{A single element numeric vector corresponding to the minimum number of simulations a contact between two residues must be present. If not, it will be set to 0 in the output matrix.
}
}
\value{
The output matrix is a nXn binary matrix (n = residue number). Elements equal to 1 correspond to residues in contact, elements equal to 0 to residues not in contact.
}
\seealso{
cmap (require bio3d package)
}
\examples{
## load example data
data(pdb)
data(trj.1)
data(trj.2)
data(trj.3)

## Lets work with Calpha atoms only
library(bio3d)
ca.inds <- atom.select(pdb, "calpha")
#noh.inds <- atom.select(pdb, "noh")

## calculate single contact map matrices
cm.1 <- cmap(trj.1[,ca.inds$xyz], pcut=0.3, scut=0, dcut=7, mask.lower=FALSE)
cm.2 <- cmap(trj.1[,ca.inds$xyz], pcut=0.3, scut=0, dcut=7, mask.lower=FALSE)
cm.3 <- cmap(trj.3[,ca.inds$xyz], pcut=0.3, scut=0, dcut=7, mask.lower=FALSE)

## create a 3D contact matrix from 3 simulations
library(abind)
cm.all <- abind(cm.1, cm.2, cm.3, along=3)

## calculate average contact matrix
cm.filter <- cmap.filter(cm=cm.all, cutoff.sims=1)

## plot the result
par(pty="s", mfcol=c(1,2))
image(cm.3, col=c(NA,"black"))
image(cm.filter, col=c(NA,"black"))
}
\keyword{analysis}
