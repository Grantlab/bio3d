\name{filter.cmap}
\alias{filter.cmap}
\title{ Contact Map Consensus Filtering }
\description{
  This function filters a tridimensional contact matrix (NxNxZ), where N is
  the residue number and Z is the simulation number) selecting only
  contacts present in at least P simulations.
}
\usage{
filter.cmap(cm, cutoff.sims = NULL)
}
\arguments{
  \item{cm}{ An array of dimensions NxNxZ or a list of NxN matrices
    containing binary contact values as obtained from
    \code{cmap}. Here, \sQuote{N} is the residue number and \sQuote{Z}
    the simulation number. The matrix elements should be 1 if two
    residues are in contact and 0 if they are not in contact. }
  \item{cutoff.sims}{ A single element numeric vector corresponding to
    the minimum number of simulations a contact between two residues
    must be present. If not, it will be set to 0 in the output matrix. }
}
\value{
  The output matrix is a nXn binary matrix (n = residue number). Elements
  equal to 1 correspond to residues in contact, elements equal to 0 to
  residues not in contact. 
}
\seealso{
  \code{\link{cmap}}, \code{\link{plot.cmap}}
}
\examples{

\dontrun{
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
  cms <- list()
  cms[[1]] <- cmap(trj1[,ca.inds$xyz], pcut=0.3, scut=0, dcut=7, mask.lower=FALSE)
  cms[[2]] <- cmap(trj1[,ca.inds$xyz], pcut=0.3, scut=0, dcut=7, mask.lower=FALSE)

  ## calculate average contact matrix
  cm.filter <- filter.cmap(cms, cutoff.sims=2)

  ## plot the result
  par(pty="s", mfcol=c(1,3))
  plot.cmap(cms[[1]])
  plot.cmap(cms[[2]])
  plot.cmap(cm.filter)
}
}
\keyword{analysis}
