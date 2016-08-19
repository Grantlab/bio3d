\name{cmap}
\alias{cmap}
\alias{cmap.default}
\alias{cmap.xyz}
\alias{cmap.pdb}
\title{ Contact Map }
\description{
  Construct a Contact Map for Given Protein Structure(s).}
\usage{
cmap(\dots)

\method{cmap}{default}(\dots)

\method{cmap}{xyz}(xyz, grpby = NULL, dcut = 4, scut = 3, pcut=1, binary=TRUE,
            mask.lower = TRUE, gc.first=FALSE, ncore=1, nseg.scale=1, \dots)

\method{cmap}{pdb}(pdb, inds = NULL, verbose = FALSE, \dots)
}
\arguments{
  \item{xyz}{ numeric vector of xyz coordinates or a numeric matrix of
    coordinates with a row per structure/frame. }
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of \code{xyz} that should be considered as a group
    (e.g. atoms from a particular residue). }  
  \item{dcut}{ a cutoff distance value below which atoms are considered
    in contact. }
  \item{scut}{ a cutoff neighbour value which has the effect of
    excluding atoms that are sequentially within this value. }
  \item{pcut}{ a cutoff probability of structures/frames showing a contact,
    above which atoms are considered in contact with respect to the ensemble.
    Ignored if \code{binary=FALSE}. }
  \item{binary}{ logical, if FALSE the raw matrix containing fraction of frames that 
    two residues are in contact is returned. }
  \item{mask.lower}{ logical, if TRUE the lower matrix elements
    (i.e. those below the diagonal) are returned as NA.}
  \item{gc.first}{ logical, if TRUE will call gc() first before calculation of 
    distance matrix. This is to solve the memory overload problem when \code{ncore > 1}
    and \code{xyz} has many rows, with a bit sacrifice on speed. }
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires package \sQuote{parallel} installed. }
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
  \item{pdb}{  a structure object of class \code{"pdb"}, obtained from
    \code{\link{read.pdb}}. }
  \item{inds}{ a list object of ATOM and XYZ indices as obtained from
    \code{\link{atom.select}}. }
  \item{verbose}{ logical, if TRUE details of the selection are
    printed. }
  \item{\dots}{ arguments passed to and from functions. }
}
\details{
  A contact map is a simplified distance matrix. See the distance matrix
  function \code{\link{dm}} for further details.
  
  Function \code{"cmap.pdb"} is a wrapper for \code{"cmap.xyz"}
  which selects all \sQuote{notwater} atoms and calculates the contact
  matrix grouped by residue number.

}
\value{
  Returns a N by N numeric matrix composed of zeros and ones, where
  one indicates a contact between selected atoms.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\seealso{ \code{\link{dm}}, \code{\link{dccm}}, \code{\link{dist}},
  \code{\link{dist.xyz}} }
\examples{

##- Read PDB file
pdb <- read.pdb( system.file("examples/hivp.pdb", package="bio3d") )

## Atom Selection indices
inds <- atom.select(pdb, "calpha")

## Reference contact map
ref.cont <- cmap( pdb$xyz[inds$xyz], dcut=6, scut=3 )
plot.cmap(ref.cont)

\dontrun{
##- Read Traj file
trj <- read.dcd( system.file("examples/hivp.dcd", package="bio3d") )
## For each frame of trajectory
sum.cont <- NULL
for(i in 1:nrow(trj)) {

  ## Contact map for frame 'i'
  cont <- cmap(trj[i,inds$xyz], dcut=6, scut=3)

  ## Product with reference
  prod.cont <- ref.cont * cont
  sum.cont <- c(sum.cont, sum(prod.cont,na.rm=TRUE))
}

plot(sum.cont, typ="l")
}
}
\keyword{ utilities }
