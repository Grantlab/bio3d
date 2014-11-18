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
  This function provides basic functionality to read and parse a AMBER
  PrmTop file. The resulting \sQuote{prmtop} object contains a complete
  list object of the information stored in the PrmTop file. 

  See examples for further details. 
}
\value{
  Returns a list of class \sQuote{prmtop} (inherits class
  \sQuote{amber}) with components according to the flags present in the PrmTop
  file. See the AMBER documentation for a complete list of
  flags/components: \url{http://ambermd.org/formats.html}.

  Selected components:
  \item{ATOM_NAME}{ a character vector of atom names. }
  \item{ATOMS_PER_MOLECULE}{ a numeric vector containing the number of
    atoms per molecule. }
  \item{MASS}{ a numeric vector of atomic masses. }
  \item{RESIDUE_LABEL}{ a character vector of residue labels. }
  \item{RESIDUE_RESIDUE_POINTER}{ a numeric vector of pointers to the
    first atom in each residue. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  \url{http://ambermd.org/formats.html}
}
\author{ Lars Skjaerven }
\note{
  See AMBER documentation for PrmTop format description:\cr
  \url{http://ambermd.org/formats.html}.
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

## Convert to multimodel PDB format
pdb <- amb2pdb(prmtop, trj[1:20,], inds=ca.inds, inds.crd=NULL)

## RMSD of trajectory
rd <- rmsd(crds$xyz[ca.inds$xyz], traj, fit=TRUE)
}
}
\keyword{ IO }
