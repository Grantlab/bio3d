\name{dm}
\alias{dm}
\alias{dm.pdb}
\alias{dm.xyz}
\alias{dm.pdbs}
\title{ Distance Matrix Analysis }
\description{
  Construct a distance matrix for a given protein structure.
}
\usage{
dm(\dots)

\method{dm}{pdb}(pdb, inds = NULL, grp = TRUE, verbose=TRUE, \dots)
\method{dm}{pdbs}(pdbs, rm.gaps=FALSE, all.atom=FALSE, 
                        aligned.atoms.only=NULL, \dots)

\method{dm}{xyz}(xyz, grpby = NULL, scut = NULL, mask.lower = TRUE, 
                  gc.first=FALSE, ncore=1, \dots)
}
\arguments{
  \item{pdb}{ a \code{pdb} structure object as returned by
    \code{\link{read.pdb}} or a numeric vector of \sQuote{xyz} coordinates.}
  \item{inds}{ atom and xyz coordinate indices obtained from \code{atom.select}
    that selects the elements of \code{pdb} upon which the calculation
    should be based. }
  \item{grp}{ logical, if TRUE atomic distances will be grouped according
    to their residue membership. See \sQuote{grpby}. }
  \item{verbose}{ logical, if TRUE possible warnings are printed. }
  \item{pdbs}{ a \sQuote{pdbs} object as returned by \code{\link{read.fasta.pdb}},
    \code{\link{read.all}}, or \code{\link{pdbaln}}. }
  \item{rm.gaps}{ logical, if TRUE gapped positions are removed in the returned value. }
  \item{all.atom}{ logical, if TRUE all-atom coordinates from \code{\link{read.all}}
    are used. }
  \item{aligned.atoms.only}{ logical, if TRUE only equivalent (aligned) atoms are considered. 
    Only meaningful when \code{all.atom=TRUE}. Default: FALSE. } 
  \item{xyz}{ a numeric vector or matrix of Cartesian coordinates.}
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of \code{xyz} that should be considered as a group
    (e.g. atoms from a particular residue). }
  \item{scut}{ a cutoff neighbour value which has the effect of
    excluding atoms, or groups, that are sequentially within this
    value.}
  \item{mask.lower}{ logical, if TRUE the lower matrix elements
    (i.e. those below the diagonal) are returned as NA.}
  \item{gc.first}{ logical, if TRUE will call gc() first before calculation of 
    distance matrix. This is to solve the memory overload problem when \code{ncore > 1}
    and \code{xyz} has many rows/columns, with a bit sacrifice on speed. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{\dots}{ arguments passed to and from functions. }
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
  distances, where N is the number of selected atoms. With multiple
  frames the output is provided in a three dimensional array. 
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
# PDB server connection required - testing excluded

##--- Distance Matrix Plot
pdb <- read.pdb( "4q21" )
k <- dm(pdb,inds="calpha")
filled.contour(k, nlevels = 10)

## NOTE: FOLLOWING EXAMPLE NEEDS MUSCLE INSTALLED
if(check.utility("muscle")) {

##--- DDM: Difference Distance Matrix
# Downlaod and align two PDB files
pdbs <- pdbaln( get.pdb( c( "4q21", "521p"), path = tempdir() ), outfile = tempfile() )

# Get distance matrix
a <- dm.xyz(pdbs$xyz[1,])
b <- dm.xyz(pdbs$xyz[2,])

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
