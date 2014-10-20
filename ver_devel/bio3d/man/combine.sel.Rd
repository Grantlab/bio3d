\name{combine.sel}
\alias{combine.sel}
\title{ Combine Atom Selections From PDB Structure }
\description{
  Do "and", "or", or "not" logical operations between two atom 
  selections made by \code{\link{atom.select}}
}
\usage{
combine.sel(sel1=NULL, sel2=NULL, op="AND", verbose=TRUE)
}
\arguments{
  \item{sel1}{ an atom selection object of class \code{"select"}, obtained from
    \code{\link{atom.select}}. }
  \item{sel2}{ a second atom selection object of class \code{"select"}, obtained 
    from \code{\link{atom.select}}. }
  \item{op}{ name of the logical operation. }
  \item{verbose}{ logical, if TRUE details of the selection combination are 
     printed. }
}
\details{
  The value of \code{op} should be one of following:
  (1) "AND", "and", or "&" for set intersection, 
  (2) "OR", "or", "|", or "+" for set union,
  (3) "NOT", "not", "!", or "-" for set difference \code{sel1 - sel2}.
}
\value{
  Returns a list of class \code{"select"} with components:
  \item{atom }{atom indices of selected atoms.}
  \item{xyz }{xyz indices of selected atoms.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao }
\seealso{ \code{\link{atom.select}} }
\examples{

# Read a PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

# Select all C-alpha atoms
ca.inds <- atom.select(pdb, "calpha")

# Select all atoms with residues numbers between 43 and 54
res.inds <- atom.select(pdb, resno=43:54)

# Select all C-beta atoms with residues numbers between 43 and 54
cb.inds <- atom.select(pdb, resno=43:54, elety="CB")

# Intersection
inds <- combine.sel(ca.inds, res.inds, op="AND")
print( pdb$atom[ inds$atom, "resid" ] )
print( pdb$xyz[ inds$xyz ] )

# Union
inds2 <- combine.sel(inds, cb.inds, op="+")

# Not
inds3 <- combine.sel(res.inds, ca.inds, op="-")

}
\keyword{utilities}

