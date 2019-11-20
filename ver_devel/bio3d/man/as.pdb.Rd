\name{as.pdb}
\alias{as.pdb}
\alias{as.pdb.mol2}
\alias{as.pdb.prmtop}
\alias{as.pdb.default}
\title{ Convert to PDB format }
\description{
  Convert Tripos Mol2 format, or Amber parameter/topology and coordinate data to PDB format.
}
\usage{
as.pdb(\dots)

\method{as.pdb}{mol2}(mol, \dots)

\method{as.pdb}{prmtop}(prmtop, crd=NULL, inds=NULL, inds.crd=inds, ncore=NULL, \dots)

\method{as.pdb}{default}(pdb=NULL, xyz=NULL, type=NULL, resno=NULL,
                   resid=NULL, eleno=NULL, elety=NULL, chain=NULL, 
                   insert=NULL, alt=NULL, o=NULL, b=NULL, segid=NULL, 
                   elesy=NULL, charge=NULL, verbose=TRUE, \dots)
}
\arguments{
  \item{\dots}{ arguments passed to and from functions. }
  \item{mol}{ a list object of type \code{"mol2"} (obtained with
    \code{\link{read.mol2}}). }
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

  \item{pdb}{ an object of class \sQuote{pdb} as obtained from
    \code{\link{read.pdb}}. }
  \item{xyz}{ a numeric vector/matrix of Cartesian coordinates. If
    provided, the number of atoms in the new PDB object will be set to
    \code{ncol(as.xyz(xyz))/3} (see \code{\link{as.xyz}}).

    If \code{xyz} is not provided the number of atoms will be based on
    the length of \code{eleno}, \code{resno}, or \code{resid} (in that
    order). }
  \item{type}{ a character vector of record types, i.e. "ATOM" or "HETATM",
    with length equal to \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element
    character vector can be provided which will be repeated to match
    the number of atoms. }
  \item{resno}{ a numeric vector of residue numbers of length equal to
    \code{ncol(as.xyz(xyz))/3}. }
  \item{resid}{ a character vector of residue types/ids of length equal to
    \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element
    character vector can be provided which will be repeated to match
    the number of atoms. }
  \item{eleno}{ a numeric vector of element/atom numbers of length equal to
    \code{ncol(as.xyz(xyz))/3}. }
  \item{elety}{ a character vector of element/atom types of length equal to
    \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element
    character vector can be provided which will be repeated to match
    the number of atoms. }
  \item{chain}{ a character vector of chain identifiers with length equal to
    \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element
    character vector can be provided which will be repeated to match
    the number of atoms. }
  \item{insert}{ a character vector of insertion code with length equal to
    \code{ncol(as.xyz(xyz))/3}. }
  \item{alt}{ a character vector of alternate record with length equal to
    \code{ncol(as.xyz(xyz))/3}. }
  \item{o}{ a numeric vector of occupancy values of length equal to
    \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element
    numeric vector can be provided which will be repeated for to match
    the number of atoms. }
  \item{b}{ a numeric vector of B-factors of length equal to \code{ncol(as.xyz(xyz))/3}.
    Alternatively, a single element numeric vector can be provided
    which will be repeated to match the number of atoms.}
  \item{segid}{ a character vector of segment id of length equal to
    \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element character vector can
    be provided which will be repeated to match the number of
    atoms. }
  \item{elesy}{ a character vector of element symbol of length equal to
    \code{ncol(as.xyz(xyz))/3}. Alternatively, a single element character vector can
    be provided which will be repeated to match the number of
    atoms. }
  \item{charge}{ a numeric vector of atomic charge of length equal to
    \code{ncol(as.xyz(xyz))/3}. }
  \item{verbose}{ logical, if TRUE details of the PDB generation process
    is printed to screen. }
}
\details{
  This function converts Tripos Mol2 format, Amber formatted parameter/topology
  (PRMTOP) and coordinate objects, and vector data to a PDB object.

  While \code{as.pdb.mol2} and \code{as.pdb.prmtop} converts specific
  objects to a PDB object, \code{as.pdb.default} provides basic
  functionality to convert raw data such as vectors of e.g. residue numbers,
  residue identifiers, Cartesian coordinates, etc to a PDB object. When
  \code{pdb} is provided the returned PDB object is built from the input
  object with fields replaced by any input vector arguments.
  e.g. \code{as.pdb(pdb, xyz=crd)} will return the same PDB object, with
  only the Cartesian coordinates changed to \code{crd}.
}
\value{
  Returns a list of class \code{"pdb"} with the following components:
  \item{atom}{ a data.frame containing all atomic coordinate ATOM data,
    with a row per ATOM and a column per record type.  See below for details
    of the record type naming convention (useful for accessing columns). }
  \item{xyz }{ a numeric matrix of ATOM coordinate data of class \code{xyz}. }
  \item{calpha }{ logical vector with length equal to \code{nrow(atom)}
    with TRUE values indicating a C-alpha \dQuote{elety}. }
  \item{call }{ the matched call. }
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
  \url{http://ambermd.org/FileFormats.php}
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{read.crd}}, \code{\link{read.ncdf}},
  \code{\link{atom.select}}, \code{\link{read.pdb}}
}
\examples{
## Vector(s) to PDB object
pdb <- as.pdb(resno=1:6, elety="CA", resid="ALA", chain="A")
pdb

\dontrun{
## Read a PRMTOP file
prmtop <- read.prmtop(system.file("examples/crambin.prmtop", package="bio3d"))

## Read Amber coordinates
crds <- read.crd(system.file("examples/crambin.inpcrd", package="bio3d"))

## Atom selection
ca.inds <- atom.select(prmtop, "calpha")

## Convert to PDB format
pdb <- as.pdb(prmtop, crds, inds=ca.inds)


## Read a single entry MOL2 file
## (returns a single object)
mol <- read.mol2( system.file("examples/aspirin.mol2", package="bio3d") )

## Convert to PDB
pdb <- as.pdb(mol)
}
}
\keyword{ IO }
