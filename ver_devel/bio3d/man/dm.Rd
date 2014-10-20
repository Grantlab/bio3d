\name{dm}
\alias{dm}
\alias{dm.xyz}
\title{ Distance Matrix Analysis }
\description{
  Construct a distance matrix for a given protein structure.
}
\usage{
dm(pdb, selection = "calpha", verbose=TRUE)
dm.xyz(xyz, grpby = NULL, scut = NULL, mask.lower = TRUE)
}
\arguments{
  \item{pdb}{ a \code{pdb} structure object as returned by
    \code{\link{read.pdb}} or a numeric vector of \sQuote{xyz} coordinates.}
  \item{selection}{ a character string for selecting the \code{pdb} atoms
    to undergo comparison (see \code{\link{atom.select}}). }
  \item{verbose}{ logical, if TRUE possible warnings are printed. }
  \item{xyz}{ a numeric vector of Cartesian coordinates.}
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of \code{xyz} that should be considered as a group
    (e.g. atoms from a particular residue). }
  \item{scut}{ a cutoff neighbour value which has the effect of
    excluding atoms, or groups, that are sequentially within this
    value.}
  \item{mask.lower}{ logical, if TRUE the lower matrix elements
    (i.e. those below the diagonal) are returned as NA.}
}
\details{
  Distance matrices, also called distance plots or distance maps, are an
  established means of describing and comparing protein conformations
  (e.g. Phillips, 1970; Holm, 1993).

  A distance matrix is a 2D representation of 3D structure that is
  independent of the coordinate reference frame and, ignoring
  chirality, contains enough information to reconstruct the 3D Cartesian
  coordinates (e.g. Havel, 1983).
}
\value{
  Returns a numeric matrix of class \code{"dmat"}, with all N by N
  distances, where N is the number of selected atoms.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.

  Phillips (1970) \emph{Biochem. Soc. Symp.} \bold{31}, 11--28.
  
  Holm (1993) \emph{J. Mol. Biol.} \bold{233}, 123--138.
  
  Havel (1983) \emph{Bull. Math. Biol.} \bold{45}, 665--720.
}
\author{ Barry Grant }
\note{
  The input \code{selection} can be any character string or pattern
  interpretable by the function \code{\link{atom.select}}.  For example,
  shortcuts \code{"calpha"}, \code{"back"}, \code{"all"} and selection
  strings of the form \code{/segment/chain/residue number/residue
    name/element number/element name/}; see \code{\link{atom.select}}
  for details.

  If a coordinate vector is provided as input (rather than a \code{pdb}
  object) the \code{selection} option is redundant and the input vector
  should be pruned instead to include only desired positions. 
}
\seealso{ \code{\link{plot.dmat}}, \code{\link{read.pdb}}, \code{\link{atom.select}} }
\examples{
\donttest{

##--- Distance Matrix Plot
pdb <- read.pdb( "4q21" )
k <- dm(pdb,selection="calpha")
filled.contour(k, nlevels = 10)

## NOTE: FOLLOWING EXAMPLE NEEDS MUSCLE INSTALLED
if(check.utility("muscle")) {

##--- DDM: Difference Distance Matrix
# Downlaod and align two PDB files
pdbs <- pdbaln( get.pdb( c( "4q21", "521p"), path = tempdir() ), outfile = tempfile() )

# Get distance matrix
a <- dm(pdbs$xyz[1,])
b <- dm(pdbs$xyz[2,])

# Calculate DDM
c <- a - b

# Plot DDM
plot(c,key=FALSE, grid=FALSE)

plot(c, axis.tick.space=10,
     resnum.1=pdbs$resno[1,],
     resnum.2=pdbs$resno[2,],
     grid.col="black",
     xlab="Residue No. (4q21)", ylab="Residue No. (521p)")
}
}

\dontrun{
##-- Residue-wise distance matrix based on the
##   minimal distance between all available atoms
l <- dm.xyz(pdb$xyz, grpby=pdb$atom[,"resno"], scut=3)
}
}
\keyword{ utilities }
