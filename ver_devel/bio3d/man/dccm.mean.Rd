\name{dccm.mean}
\alias{dccm.mean}
\title{
Atomic Correlation Matrix Consensus Filtering
}
\description{
This function filters a tridimensional correlation matrix (nXnXz, where n is the residue number and z is the simulation number) selecting only the correlations present in at least p simulations, where p<=z.
}
\usage{
dccm.mean(cij, cutoff.sims = dim(cij)[3], cutoff.cij = 0.4)
}
\arguments{
  \item{cij}{ A numeric array with 3 dimensions (nXnXz) containing atomic correlation values (cij), where 
"n" is the residue number, "z" the simulation number. cij values must be between -1 and 1 (see Pearson's correlation).  
}
  \item{cutoff.sims}{ A single element numeric vector corresponding to the minimum number of simulations a correlation between two residues must be present and whose absolute value must be higher than cutoff.cij. If not, the matrix element will be set to 0 and not used for the average calculation.
}
  \item{cutoff.cij}{  A single element numeric vector corresponding to the minimum absolute value a correlation element must be equal to or higher to be used for the average calculation. If lower, the element will be set to 0 and not used for the average calculation.
}
}
\value{
This function returns a numeric atomic correlation matrix. Classes: dccm, matrix
}

\seealso{
dccm (require bio3d package)
}
\examples{
###############################
######## Play example #########
###############################

pth <- system.file("data",package="cna")  ## path to the example files
load(paste(pth,"/play.cij.RData",sep=""))

library(bio3d)
library(abind)

cij.all <- abind(cij.play.1, cij.play.2, cij.play.3, along=3)

dccm.mean(cij.all,cutoff.sims=1)
dccm.mean(cij.all,cutoff.sims=3)

dccm.mean(cij.all,cutoff.sims=2,cutoff.cij=0.6)
dccm.mean(cij.all,cutoff.sims=2,cutoff.cij=0.4)


###############################
## load protein example data ##
###############################

data(pdb)
data(trj.1)
data(trj.2)
data(trj.3)

## Lets work with Calpha atoms only
library(bio3d)
ca.inds <- atom.select(pdb, "calpha")
#noh.inds <- atom.select(pdb, "noh")

## generation of the cij for each simulations (3 sims)
cij.1 <- dccm( trj.1[,ca.inds$xyz] )
cij.2 <- dccm( trj.2[,ca.inds$xyz] )
cij.3 <- dccm( trj.3[,ca.inds$xyz] )

## generation of the 3D cij matrix
library(abind)
cij.all <- abind(cij.1, cij.2, cij.3, along=3)

## calculation of the average cij matrix. To be used in the average calculation
## each contact must have an absolute value higher than 0.4 in at least 2 sims.
cij.ave <- dccm.mean(cij.all,2)

## plot the result
par(pty="s", mfcol=c(1,2))
image(cij.1, col=c(NA,"black"))    ## For simplicity cij values are plotted all black
image(cij.ave, col=c(NA,"black"))  ## For simplicity cij values are plotted all black

}
\keyword{analysis}

