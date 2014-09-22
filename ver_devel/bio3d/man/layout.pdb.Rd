\name{layout.pdb}
\alias{layout.pdb}
\title{
 Geometric Center Calculator
}
\description{
  Determine the 3D geometric center of each community
  from its C-alpha atom positions.
}
\usage{
  layout.pdb(pdb, membership, renumber=TRUE, k=3)
}
\arguments{
  \item{pdb}{ A pdb class object as obtained from \sQuote{read.pdb} function. }
  \item{membership}{ A numeric vector with community residue membership or a 
  network object as obtained from {cna} function. }
  \item{renumber}{ Logical, if TRUE the \sQuote{pdb} will be renumberd starting
  at resisidue number one before community coordinate averages are calculated. 
  Care needs to be taken to ensure your residue numbers in the input pdb object 
  and comunity vector names match. }
  \item{k}{ the maximum dimension of the space which the data are to be
          represented in; must be between 1 and 3. }
}
\details{
  This function calculates the geometric center for each community from
  the atomic position of it's Calpha atoms taken from a pdb file.
  The community residue membership are given as input.
}
\value{
 A numeric matrix of Nxk, where N is the number of communities and k the number of dimensions requested.
}
\author{ Guido Scarabelli and Barry Grant }

\seealso{ \code{\link{plot.igraph}}, \code{\link{plot.communities}},
  \code{\link{igraph.plotting}} \code{\link{plot.cna}}}
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

# Plot will be slow
xy <- plot.cna(net)
plot3d.cna(net, pdb)

layout.pdb(pdb, net, k=3)
layout.pdb(pdb, net, k=2)

# can be used as input to plot.cna and plot3d.cna....
# plot.cna(net, layout=layout.pdb(pdb, net, k=2))
# plot3d.cna(net, pdb, layout=layout.pdb(pdb, net, k=2))

}
\keyword{ utility }