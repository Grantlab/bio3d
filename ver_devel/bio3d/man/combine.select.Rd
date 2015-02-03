\name{combine.select}
\alias{combine.select}
\title{ Combine Atom Selections From PDB Structure }
\description{
  Do "and", "or", or "not" set operations between two or more atom
  selections made by \code{\link{atom.select}}
}
\usage{
combine.sel(sel1=NULL, sel2=NULL, \dots, operator="AND", verbose=TRUE)
}
\arguments{
  \item{sel1}{ an atom selection object of class \code{"select"}, obtained from
    \code{\link{atom.select}}. }
  \item{sel2}{ a second atom selection object of class \code{"select"}, obtained
    from \code{\link{atom.select}}. }
  \item{\dots}{ more select objects for the set operation. }
  \item{operator}{ name of the set operation. }
  \item{verbose}{ logical, if TRUE details of the selection combination are
     printed. }
}
\details{
  The value of \code{operator} should be one of following:
  (1) "AND", "and", or "&" for set intersect,
  (2) "OR", "or", "|", or "+" for set union,
  (3) "NOT", "not", "!", or "-" for set difference \code{sel1 - sel2 - sel3 ...}.
}
\value{
  Returns a list of class \code{"select"} with components:
  \item{atom }{atom indices of selected atoms.}
  \item{xyz }{xyz indices of selected atoms.}
  \item{call }{the matched call.}
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao }
\seealso{ \code{\link{atom.select}} }
\examples{

# Read a PDB file
pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )

## - Build atom selections to be operated
# Select C-alpha atoms of entire system
ca.global.inds <- atom.select(pdb, "calpha")

# Select C-beta atoms of entire protein
cb.global.inds <- atom.select(pdb, "protein", elety="CB")

# Select backbone atoms of entire system
bb.global.inds <- atom.select(pdb, "backbone")

# Select all atoms with residue number from 46 to 50
aa.local.inds <- atom.select(pdb, resno=46:50)


# Do set intersect: return C-alpha atoms with residue number from 46 to 50
ca.local.inds <- combine.sel(ca.global.inds, aa.local.inds)
print( pdb$atom[ ca.local.inds$atom, ] )

# Do set subtract: return side-chain atoms with residue number from 46 to 50
sc.local.inds <- combine.sel(aa.local.inds, bb.global.inds, op="-")
print( pdb$atom[ sc.local.inds$atom, ] )

# Do set union: return C-alpha and side-chain atoms with residue number from 46 to 50
casc.local.inds <- combine.sel(ca.local.inds, sc.local.inds, op="+")
print( pdb$atom[ casc.local.inds$atom, ] )

# More than two selections: return side-chain atoms (but not C-beta) with residue number from 46 to 50
sc2.local.inds <- combine.sel(aa.local.inds, bb.global.inds, cb.global.inds, op="-")
print( pdb$atom[ sc2.local.inds$atom, ] )

}
\keyword{utilities}

