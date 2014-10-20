\encoding{UTF-8}
\name{atom2mass}
\alias{atom2mass}
\alias{atom2mass.default}
\alias{atom2mass.pdb}
\title{ Atom Names/Types to Mass Converter }
\description{
  Convert atom names/types into atomic masses.
}
\usage{
atom2mass(\dots)
\method{atom2mass}{default}(x, mass.custom=NULL, elety.custom=NULL,
                            grpby=NULL, rescue=TRUE, \dots) 
\method{atom2mass}{pdb}(pdb, inds=NULL, mass.custom=NULL,
                        elety.custom=NULL, grpby=NULL, rescue=TRUE, \dots) 
}
\arguments{
  \item{x}{ a character vector containing atom names/types to be
    converted. } 
  \item{mass.custom}{ a customized data.frame containing atomic symbols
    and corresponding masses. } 
  \item{elety.custom}{ a customized data.frame containing atom
    names/types and corresponding atomic symbols.} 
  \item{grpby}{a \sQuote{factor}, as returned by \code{as.factor}, used
    to group the atoms.} 
  \item{rescue}{ logical, if TRUE the atomic symbols will be mapped to
    the first character of the atom names/types.} 
  \item{pdb}{ an object of class \sQuote{pdb} for which \code{elety} will
    be converted.} 
  \item{inds}{ an object of class \sQuote{select} indicating a subset of
    the \code{pdb} object to be used (see \code{\link{atom.select}} and
    \code{\link{trim.pdb}}).} 
  \item{\dots}{.}
}
\details{
  The default method first convert atom names/types into atomic symbols
  using the \code{\link{atom2ele}} function. Then, atomic symbols are
  searched in the \code{elements} data set and their corresponding masses
  are returned. If \code{mass.custom} is specified it is combined with
  \code{elements} (using \code{rbind}) before searching. Therefore,
  \code{mass.custom} must have columns named \code{symb} and \code{mass}
  (see examples). If \code{grpby} is specified masses are splitted (using
  \code{split}) to compute the mass of groups of atoms defined by
  \code{grpby}. 

  The S3 method for object of class \sQuote{pdb}, pass
  \code{pdb$atom$elety} to the default method. 
}
\value{Return a numeric vector of masses.}
\author{Julien Ide, Lars Skjaerven}
\seealso{
  \code{\link{elements}}, \code{\link{atom.index}},
  \code{\link{atom2ele}}, \code{\link{read.pdb}}
}
\examples{
atom.names <- c("CA", "O", "N", "OXT")
atom2mass(atom.names)

\donttest{

## Get atomic symbols from a PDB object with a customized data set
pdb <- read.pdb("3RE0", verbose=FALSE)
inds <- atom.select(pdb, resno=201, verbose=FALSE)

## selected atoms
print(pdb$atom$elety[inds$atom])

## default will map CL2 to C
atom2mass(pdb, inds)

## map element CL2 correctly to Cl
myelety  <- data.frame(name = c("CL2","PT1","N1","N2"), symb = c("Cl","Pt","N","N"))
atom2mass(pdb, inds, elety.custom = myelety)

## custom masses
mymasses <- data.frame(symb = c("Cl","Pt"), mass = c(35.45, 195.08))
atom2mass(pdb, inds, elety.custom = myelety, mass.custom = mymasses)

}
}
\keyword{ utilities }
