\name{mask}
\alias{mask}
\alias{mask.dccm}
\title{ Mask a Subset of Atoms in a DCCM Object. }
\description{
  Produce a new DCCM object with selected atoms masked.
}
\usage{
mask(\dots)

\method{mask}{dccm}(dccm, pdb = NULL, a.inds = NULL, b.inds = NULL, \dots)
}
\arguments{
  \item{dccm}{ a DCCM structure object obtained from function
    \code{\link{dccm}}. }
  \item{pdb}{ a PDB structure object obtained from
    \code{\link{read.pdb}}. Must match the dimensions of \code{dccm}. }
  \item{a.inds}{ a list object of ATOM and XYZ indices as obtained from
    \code{\link{atom.select}} specifying which elements of the DCCM
    matrix which should not be masked. }
  \item{b.inds}{ a list object of ATOM and XYZ indices as obtained from
    \code{\link{atom.select}} specifying which elements of the DCCM
    matrix which should not be masked. }
  \item{\dots}{ arguments not passed anywhere. }
}
\details{
  This is a basic utility function for masking a DCCM object matrix to
  highlight a few user-selected areas in the correlation network.

  Note: The current version assumes that the input PDB corresponds to
  the input DCCM. In many cases this will correspond to a PDB
  object containing only CA atoms.
}
\value{
  Returns a matrix list of class \code{"dccm"} with the indices/atoms
  not corresponding to the selection masked.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{dccm}}, \code{\link{atom.select}}
}
\examples{
## Calculate DCCM
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
cij <- dccm(nma(pdb))

## Mask DCCM matrix according to matrix indices
cijm <- mask(cij, a.inds=40:50, b.inds=80:90)
plot(cijm)


## Mask DCCM matrix according PDB selection
pdb.ca <- trim(pdb, "calpha")
a.inds <- atom.select(pdb.ca, resno=40:50)
b.inds <- atom.select(pdb.ca, resno=80:90)

# Provide pdb object correspoding to input dccm
cijm <- mask(cij, pdb.ca, a.inds, b.inds)
plot(cijm)
}
\keyword{ utilities }
