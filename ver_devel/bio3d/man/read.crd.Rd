\name{read.crd}
\alias{read.crd}
\title{ Read Coordinate Data from Amber or Charmm }
\description{
  Read a CHARMM CARD (CRD) or AMBER coordinate file.
}
\usage{
read.crd(file, ...)
}
\arguments{
  \item{file}{ the name of the coordinate file to be read.  }
  \item{\dots}{ additional arguments passed to the methods
    \code{read.crd.charmm} or \code{read.crd.amber}. }
}
\details{
  \code{read.crd} is a generic function calling the corresponding function
  determined by the class of the input argument \code{x}. Use
  \code{methods("read.crd")} to get all the methods for \code{read.crd}
  generic:

  \code{\link{read.crd.charmm}} will be used for file extension
  \sQuote{.crd}.

  \code{\link{read.crd.amber}} will be used for file extension
  \sQuote{.rst} or \sQuote{.inpcrd}.

  See examples for each corresponding function for more details.
}
\value{
  See the \sQuote{value} section for the corresponding functions for
  more details.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant and Lars Skjaerven}
\seealso{
  \code{\link{read.crd.amber}}, \code{\link{read.crd.charmm}},
  \code{\link{write.crd}}, \code{\link{read.prmtop}},
  \code{\link{read.pdb}}, \code{\link{write.pdb}},
  \code{\link{atom.select}},
  \code{\link{read.dcd}}, \code{\link{read.ncdf}} }
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
