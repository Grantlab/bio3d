\name{atom.select.prmtop}
\alias{atom.select.prmtop}
\title{ Atom Selection From a Amber PrmTop object }
\description{
  Return the atom and xyz coordinates indices of a \sQuote{prmtop} structure
  object corresponding to the intersection of a hierarchical selection.
}
\usage{
\method{atom.select}{prmtop}(prmtop, ...)
}
\arguments{
  \item{prmtop}{ a structure object of class \code{"prmtop"}, obtained from
    \code{\link{read.prmtop}}. }
  \item{\dots}{ arguments passed to method \code{\link{atom.select.pdb}}. } 
}
\details{
  See \code{\link{atom.select.pdb}} for details.
}
\value{
  Returns a list of class \code{"select"} with components:
  \item{atom }{atom indices of selected atoms.}
  \item{xyz }{xyz indices of selected atoms.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Lars Skjaerven }
\seealso{
  \code{\link{atom.select.pdb}}, \code{\link{combine.sel}},
  \code{\link{read.prmtop}} 
}
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
\keyword{utilities}

