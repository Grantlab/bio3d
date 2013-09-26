\name{read.dcd}
\alias{read.dcd}
\title{ Read CHARMM/X-PLOR/NAMD Binary DCD files }
\description{
  Read coordinate data from a binary DCD trajectory file.
}
\usage{
read.dcd(trjfile, big=FALSE, verbose = TRUE)
}
\arguments{
  \item{trjfile}{ name of trajectory file to read. A vector if 
    treat a batch of files }
  \item{big}{ logical, if TRUE attempt to read large files into a
    big.matrix object }
  \item{verbose}{ logical, if TRUE print details of the reading
    process. }
}
\details{
  Reads a CHARMM or X-PLOR/NAMD binary trajectory file with either big-
  or little-endian storage formats.
  
  Reading is accomplished with two different sub-functions:
  \code{dcd.header}, which reads header info, and \code{dcd.frame}, which
  takes header information and reads atoms frame by frame producing an
  nframes/natom*3 matrix of cartesian coordinates.

}
\value{
   A numeric matrix of xyz coordinates with a frame/structure per row and
   a Cartesian coordinate per column.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant }
\note{
  See CHARMM documentation for DCD format description.

  If you experience problems reading your trajectory file with read.dcd()
  consider first reading your file into VMD and from there exporting a
  new DCD trajectory file with the 'save coordinates' option. This new
  file should be easily read with read.dcd().

   Error messages beginning 'cannot allocate vector of size' indicate a
   failure to obtain memory, either because the size exceeded the
   address-space limit for a process or, more likely, because the system
   was unable to provide the memory. Note that on a 32-bit OS there may
   well be enough free memory available, but not a large enough
   contiguous block of address space into which to map it. In such cases
   try setting the input option 'big' to TRUE. This is an experimental
   option that results in a 'big.matrix' object.

}
\seealso{ \code{\link{read.pdb}}, \code{\link{write.pdb}},
  \code{\link{atom.select}} }
\examples{
##-- Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

## Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb,"///24:27,85:90///CA/")

## lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

##-- RMSD of trj frames from PDB
r1 <- rmsd(a=pdb, b=xyz)

\dontrun{
# Pairwise RMSD of trj frames for positions 47 to 54
flap.inds <- atom.select(pdb,"///47:54///CA/")
p <- rmsd(xyz[,flap.inds$xyz])
# plot highlighting flap opening?
plot.dmat(p, color.palette = mono.colors)
}

}
\keyword{ IO }
