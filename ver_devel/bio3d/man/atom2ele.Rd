\name{atom2ele}
\alias{atom2ele}
\alias{atom2ele.default}
\alias{atom2ele.pdb}
\title{ Atom Names/Types to Atomic Symbols Converter }
\description{
  Convert atom names/types into atomic symbols
}
\usage{
atom2ele(\dots)

\method{atom2ele}{default}(x, elety.custom=NULL, rescue=TRUE, \dots)

\method{atom2ele}{pdb}(pdb, inds=NULL, elety.custom=NULL, rescue=TRUE, \dots)
}
\arguments{
  \item{x}{a character vector containing atom names/types to be converted.}
  \item{elety.custom}{a customized data.frame containing atom
    names/types and corresponding atomic symbols.} 
  \item{rescue}{logical, if TRUE the atomic symbols will be converted 
    based on matching with \code{bio3d::elements$symb}. }
  \item{pdb}{an object of class \sQuote{pdb} for which \code{elety} will
    be converted.} 
  \item{inds}{an object of class \sQuote{select} indicating a subset of
    the \code{pdb} object to be used (see \code{\link{atom.select}} and
    \code{\link{trim.pdb}}).} 
  \item{\dots}{further arguments passed to or from other methods.} 
}
\details{
  The default method searchs for the atom names/types in the
  \code{\link{atom.index}} data set and returns their corresponding atomic
  symbols. If \code{elety.custom} is specified it is combined with
  \code{\link{atom.index}} (using \code{rbind}) before
  searching. Therefore, \code{elety.custom} must contains columns named
  \code{name} and \code{symb}. 
  
  The S3 method for object of class \sQuote{pdb}, pass
  \code{pdb$atom[,"elety"]} to the default method. 
}
\value{Return a character vector of atomic symbols}
\author{Julien Ide, Lars Skjaerven}
\seealso{
  \code{\link{atom.index}},  \code{\link{elements}},
  \code{\link{read.pdb}},
  \code{\link{atom2mass}}, \code{\link{formula2mass}}
}
\examples{
atom.names <- c("CA", "O", "N", "OXT")
atom2ele(atom.names)

\donttest{
# PDB server connection required - testing excluded

## Get atomic symbols from a PDB object with a customized data set
pdb <- read.pdb("3RE0",verbose=FALSE)
lig <- trim(pdb, "ligand")

## maps PT1 to Pt, CL2 to Cl, C4A to C
atom2ele(lig)

## map atom name to element manually
myelety <- data.frame(name = "CL2", symb = "Cl")
atom2ele(lig, elety.custom = myelety)

}
}
\keyword{ utilities }
