\name{filter.cmap}
\alias{filter.cmap}
\title{ Contact Map Consensus Filtering }
\description{
This function filters a tridimensional contact matrix (nXnXz, where n is the residue number and z is the simulation number) selecting only the contact present in at least p simulations, where p<=z.
}
\usage{
filter.cmap(cm, cutoff.sims = dim(cm)[3])
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
  \code{\link{cmap}}
}
\examples{

\dontrun{
  ## need abind package
  if(!require(abind)) {
     install.packages("abind")
     require(abind)
  }

  ## load example data
  pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
  pdb <- read.pdb(pdbfile)

  trtfile <- system.file("examples/hivp.dcd", package="bio3d")
  trj <- read.dcd(trtfile, verbose=FALSE)

  ## split the trj example in two
  num.of.frames <- dim(trj)[1]
  trj1 <- trj[1:(num.of.frames/2),]
  trj2 <- trj[((num.of.frames/2)+1):num.of.frames,]

  ## Lets work with Calpha atoms only
  ca.inds <- atom.select(pdb, "calpha")
  #noh.inds <- atom.select(pdb, "noh")

  ## calculate single contact map matrices
  cm.1 <- cmap(trj1[,ca.inds$xyz], pcut=0.3, scut=0, dcut=7, mask.lower=FALSE)
  cm.2 <- cmap(trj2[,ca.inds$xyz], pcut=0.3, scut=0, dcut=5, mask.lower=FALSE)

  ## create a 3D contact matrix from 3 simulations
  cm.all <- abind(cm.1, cm.2, along=3)

  ## calculate average contact matrix
  cm.filter <- filter.cmap(cm=cm.all, cutoff.sims=2)

  ## plot the result
  par(pty="s", mfcol=c(1,3))
  image(cm.1, col=c(NA,"black"))
  image(cm.2, col=c(NA,"black"))
  image(cm.filter, col=c(NA,"black"))
}
}
\keyword{analysis}
