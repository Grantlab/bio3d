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
  \item{\dots}{ arguments passed to and from functions. }
}
\details{
  Read a AMBER Coordinate format file.
}
\value{
  A list object of type \sQuote{amber} and \sQuote{crd} with the
  following components:
  \item{xyz}{ a numeric matrix of class \sQuote{xyz} containing the
    Cartesian coordinates. }
  \item{velocities}{ a numeric vector containg the atom velocities. }
  \item{time}{ numeric, length of the simulation (applies to Amber
    restart coordinate files). }
  \item{natoms}{ total number of atoms in the coordinate file. }
  \item{box}{ dimensions of the box. }
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
  \code{\link{as.pdb}}, \code{\link{atom.select}},
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
pdb <- as.pdb(prmtop, crds, inds=ca.inds)
}
}
\keyword{ IO }
