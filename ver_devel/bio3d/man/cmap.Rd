\name{cmap}
\alias{cmap}
\title{ Contact Map }
\description{
  Construct a Contact Map for Given Protein Structure(s).}
\usage{
cmap(xyz, grpby=NULL, dcut = 4, scut = 3, pcut=1, mask.lower = TRUE,
         ncore=1, nseg.scale=1)
}
\arguments{
  \item{xyz}{ numeric vector of xyz coordinates or a numeric matrix of
    coordinates with a row per structure/frame. }
  \item{grpby}{ a vector counting connective duplicated elements that
    indicate the elements of \code{xyz} that should be considered as a group
    (e.g. atoms from a particular residue). }  
  \item{dcut}{ a cuttoff distance value below which atoms are considered
    in contact. }
  \item{scut}{ a cuttoff neighbour value which has the effect of
    excluding atoms that are sequentially within this value. }
  \item{pcut}{ a cutoff percentage of structures/frames showing a contact,
    above which atoms are considered in contact with respect to the ensemble}
  \item{mask.lower}{ logical, if TRUE the lower matrix elements
    (i.e. those below the diagonal) are returned as NA.}
  \item{ncore }{ number of CPU cores used to do the calculation.
    \code{ncore>1} requires multicore package installed. }
  \item{nseg.scale }{ split input data into specified number of segments
    prior to running multiple core calculation. See \code{\link{fit.xyz}}. }
}
\details{
  A contact map is a simplified distance matrix. See the distance matrix
  function \code{\link{dm}} for further details.
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

##- Read PDB and Traj files
pdb <- read.pdb( system.file("examples/hivp.pdb", package="bio3d") )
##trj <- read.dcd( system.file("examples/hivp.dcd", package="bio3d") )

## Atom Selection indices
inds <- atom.select(pdb, "calpha")

## Reference contact map
ref.cont <- cmap( pdb$xyz[inds$xyz], dcut=6, scut=3 )
plot.dmat(ref.cont)

\dontrun{
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
