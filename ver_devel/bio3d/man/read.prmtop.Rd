\name{read.prmtop}
\alias{read.prmtop}
\title{ Read AMBER Parameter/Topology files }
\description{
  Read parameter and topology data from an AMBER PrmTop file.
}
\usage{
read.prmtop(file)
}
\arguments{
  \item{file}{ name of prmtop file to read. }
}
\details{
  Reads a AMBER PrmTop format file.
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
  See AMBER documentation for PrmTop format description.
}
\seealso{
  \code{\link{read.crd}}, \code{\link{read.ncdf}},
  \code{\link{amb2pdb}}, \code{\link{atom.select}},
  \code{\link{read.pdb}} }
\examples{
\dontrun{
## Read a PRMTOP file
prmtop <- read.prmtop("prot_solvated.prmtop")
print(prmtop)

## Explore parameters
prmtop$MASS
prmtop$ATOM_NAME

## Atom selection 
ca.inds <- atom.select(prmtop, "calpha")

## Trajectory processing
trj <- read.ncdf("traj.nc", at.sel=ca.inds)

## Read a Amber CRD file
crds <- read.crd("prot_solvated.inpcrd")

## Convert to PDB format
pdb <- amb2pdb(prmtop, crds, inds=ca.inds)

## Convert to PDB format
pdb <- amb2pdb(prmtop, trj[20,], inds=ca.inds, inds.crd=NULL)


}
}
\keyword{ IO }
