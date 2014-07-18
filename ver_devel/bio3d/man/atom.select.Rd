\name{atom.select}
\alias{atom.select}
\title{ Atom Selection }
\description{
  Return the atom and xyz coordinates indices of a \sQuote{pdb} or Amber
  \sQuote{prmtop} structure object corresponding to the intersection of
  a hierarchical selection. 
}
\usage{
atom.select(...)
}
\arguments{
  \item{\dots}{ arguments passed to methods
    \code{\link{atom.select.pdb}} or \code{\link{atom.select.prmtop}}. } 
}
\details{
  \code{atom.select} is a generic function calling the corresponding function
  determined by the class of the input argument \code{x}. Use
  \code{methods("atom.select")} to get all the methods for \code{atom.select}
  generic:
  
  \code{\link{atom.select.pdb}} will be used when \code{x} is an object of 
  \code{"pdb"}.
  
  \code{\link{atom.select.prmtop}} will be used when \code{x} is an object of 
  \code{"prmtop"}.
  
  See examples for each corresponding function for more details.
}
\value{
  Returns a list of class \code{"select"} with components:
  \item{atom }{atom indices of selected atoms.}
  \item{xyz }{xyz indices of selected atoms.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ 
  \code{\link{combine.sel}}, \code{\link{read.pdb}}, 
  \code{\link{write.pdb}}, \code{\link{read.dcd}} 
}
\examples{

# Read a PDB file
pdb <- read.pdb( "http://www.rcsb.org/pdb/files/1BG2.pdb" )

# Print structure summary
atom.select(pdb)

# Select all C-alpha atoms with residues numbers between 65 and 143
ca.inds <- atom.select(pdb, resno=65:143, elety="CA")
print( pdb$atom[ ca.inds$atom, "resid" ] )
print( pdb$xyz[ ca.inds$xyz ] )

\dontrun{
# Select all C-alphas
ca.inds <- atom.select(pdb, "calpha")

# String examples (see above).
ca.inds <- atom.select(pdb, "///65:143///CA/")
inds <- atom.select(pdb, "///130:142///N,CA,C,O/")
}
}
\keyword{utilities}

