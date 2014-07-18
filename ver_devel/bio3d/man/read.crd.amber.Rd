\name{read.crd.amber}
\alias{read.crd.amber}
\title{ Read AMBER Coordinate files }
\description{
  Read coordinate data from an AMBER coordinate / restart file.
}
\usage{
\method{read.crd}{amber}(file, ...)
}
\arguments{
  \item{file}{ name of crd file to read. }
  \item{\dots}{ arguments going nowhere. } 
}
\details{
  Reads a AMBER Coordinate format file.
}
\value{
  A list of stuff... 
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  \url{http://ambermd.org/formats.html}
}
\author{ Lars Skjaerven }
\note{
  See AMBER documentation for Coordinate format description.
}
\seealso{
  \code{\link{read.prmtop}}, \code{\link{read.ncdf}},
  \code{\link{amb2pdb}}, \code{\link{atom.select}},
  \code{\link{read.pdb}}, \code{\link{read.crd.charmm}}  }
\examples{
\dontrun{
## Read a PRMTOP file
prmtop <- read.prmtop("prot_solvated.prmtop")
print(prmtop)

## Read a Amber CRD file
crds <- read.crd("prot_solvated.inpcrd")

## Atom selection 
ca.inds <- atom.select(prmtop, "calpha")

## Convert to PDB format
pdb <- amb2pdb(prmtop, crds, inds=ca.inds)
}
}
\keyword{ IO }
