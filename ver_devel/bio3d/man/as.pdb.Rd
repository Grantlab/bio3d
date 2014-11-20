\name{as.pdb}
\alias{as.pdb}
\alias{as.pdb.prmtop}
\title{ Amber PRMTOP to PDB converter }
\description{
  Convert Amber parameter/topology and coordinate data to PDB format.
}
\usage{
as.pdb(\dots)

\method{as.pdb}{prmtop}(pdb, crd, inds = NULL, inds.crd = inds, ncore=NULL, \dots)

}
\arguments{
  \item{\dots}{ arguments passed to and from functions. }
  \item{prmtop}{ a list object of type \code{"prmtop"} (obtained with
    \code{\link{read.prmtop}}). }
  \item{crd}{ a list object of type \code{"crd"} (obtained with
    \code{\link{read.crd.amber}}). }
  \item{inds}{ a list object of type \code{"select"} as obtained from
    \code{\link{atom.select}}. The indices points to which atoms in the
    PRMTOP object to convert. }
  \item{inds.crd}{ same as the \sQuote{inds} argument, but pointing to
    the atoms in CRD object to convert. By default, this argument equals
    to \sQuote{inds}, assuming the same number and sequence of atoms in the
    PRMTOP and CRD objects. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
}
\details{
  Converts Amber formatted parameter/topology (PRMTOP) and coordinate
  objects to a PDB object. 
}
\value{
  Returns a list of class \code{"pdb"} with the following components:
  \item{atom}{ a character matrix containing all atomic coordinate ATOM data, 
    with a row per ATOM and a column per record type.  See below for details
    of the record type naming convention (useful for accessing columns). }
  \item{xyz }{ a numeric matrix of ATOM coordinate data. }
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  \url{http://ambermd.org/formats.html}
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{read.crd}}, \code{\link{read.ncdf}},
  \code{\link{atom.select}}, \code{\link{read.pdb}}
}
\examples{
\dontrun{
## Read a PRMTOP file
prmtop <- read.prmtop("prot_solvated.prmtop")

## Atom selection 
ca.inds <- atom.select(prmtop, "calpha")

## Read a Amber CRD file
crds <- read.crd("prot_solvated.inpcrd")

## Convert to PDB format
pdb <- as.pdb(prmtop, crds, inds=ca.inds)

## Convert to PDB format
trj <- read.ncdf("traj.nc", at.sel=ca.inds)
pdb <- as.pdb(prmtop, trj[20,], inds=ca.inds, inds.crd=NULL)
}
}
\keyword{ IO }
