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

##--- Distance Matrix Plot
pdb <- read.pdb( system.file("examples/d1bg2__.ent", package = "bio3d") )
k <- dm(pdb,selection="calpha")
filled.contour(k, nlevels = 4)


##--- DDM: Difference Distance Matrix
# Read aligned PDBs
aln <- read.fasta(system.file("examples/kif1a.fa",package="bio3d"))
pdb.path=paste(system.file(package="bio3d"),"/examples/",sep="")
m <- read.fasta.pdb(aln, pdb.path = pdb.path, pdbext = ".ent")

# Get distance matrix
a <- dm(m$xyz[2,])
b <- dm(m$xyz[3,])

# Calculate DDM
c <- a - b

# Plot DDM
plot(c,key=FALSE, grid=FALSE)

plot(c, axis.tick.space=10,
     resnum.1=m$resno[1,],
     resnum.2=m$resno[2,],
     grid.col="gray",
     xlab="Residue No. (1i6i)", ylab="Residue No. (1i5s)")

\dontrun{
##-- Residue-wise distance matrix based on the
##   minimal distance between all available atoms
l <- dm.xyz(pdb$xyz, grpby=pdb$atom[,"resno"], scut=3)


##--- Extract all-atom contacts
pdb <- read.pdb( system.file("examples/d1bg2__.ent", package = "bio3d") )
l <- dm(pdb,selection="all")
l[upper.tri(l)]=NA  # make top diagonal NA

# Find residues with contacting atoms (<=5 Angstrom)
inds.stru <- which(l<=5, arr.ind=TRUE)   

# Find non-consecutive residues (>5 residues sequence separation)
seq.sep  <- abs(as.numeric(pdb$atom[inds.stru[,1],"resno"]) -
               as.numeric(pdb$atom[inds.stru[,2],"resno"]))
inds.seq <- which(seq.sep>5) # seperated by > 5 residues

# All-atom contacts (indices)
inds <- inds.stru[inds.seq,]

# All-atom contacts (now in terms of residue numbers)
tmp <- unique( paste(pdb$atom[inds[,1],"resno"],
                     pdb$atom[inds[,2],"resno"], sep="#") )

contacts <- matrix(as.numeric(unlist(strsplit(tmp,split="#"))),
                   ncol=2, byrow=TRUE )

# Plot residue contacts
freq  <- table(contacts)
xaxis <-as.vector(bounds(as.numeric(names(freq)))[,c(1,2)])
x11()
plot(freq, typ="h", xlab="Residue Number",
     xaxt="n",ylab="Number of Contacts" )
axis(1,at=xaxis,labels=xaxis)
}
}
\keyword{ utilities }
