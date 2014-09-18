\name{prune.couplings}
\alias{prune.couplings}
\title{
  Removal of couplings among residues belonging to the same
  secondary structure element.
}
\description{
  This function removes the couplings (i.e. network edges) among X
  neighboring residues forming an alpha helix or a beta strand. It is
  useful to fast WISP suboptimal paths calculation.
}
\usage{
  prune.couplings(x, sse, neigh=2)
}
\arguments{
  \item{x}{ A numeric matrix used to build a \sQuote{cna} network. }
  \item{sse}{ A secondary structure object such as obtained with
    \sQuote{dssp} function. }
  \item{neigh}{ Numeric. Specify the number of neighboring residues \bold{from}
    which the couplings are removed (i.e. neigh=2 removes couplings
    between res1 and res3, res1 and res4, etc, not touching the value
    among res1 and res2).}
}
\value{
  Returns a numeric matrix with the pruned couplings that can be used as
  input with WISP.
}
\details{
  The input of this function must be an adjacency matrix used to build a
  network with the \sQuote{cna} function. Removing the edges inside each
  secondary structural elements speeds up a bit the calculation of the
  suboptimal paths with WISP program.
}  
\seealso{
   \code{\link{cna}}, WISP: http://nbcr.ucsd.edu/data/sw/hosted/wisp
}
\examples{

# Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb,"///24:27,85:90///CA/")

# lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

# calculate dynamical cross-correlation matrix
cij <- dccm(xyz)
#plot(cij)

# Build, and betweeness cluster, a network graph
net <- cna(cij)

\dontrun{
## the example is missing sse
sse <- dssp(sse)
pruned.matrix <- prune.couplings(net$cij, sse=sse, neigh=2)

write.table(pruned.matrix, quote=FALSE, row.names=FALSE,
col.names=FALSE, file="pruned.matrix.txt")
}
}
\keyword{utilities}
