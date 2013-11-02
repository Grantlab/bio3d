\name{membership.comparison}
\alias{membership.comparison}
\title{
  Community Comparison Plot 
}
\description{
  This function compares the community clustering of a specific subset
  of residues obtained with different parameters.
}
\usage{
  membership.comparison(x, important.members=important.members, plot=FALSE, comm.color=1, sse=NULL, dash.lines=FALSE)
}
\arguments{
  \item{x}{ A numeric matrix with MxN dimensions, where M is equal to
    the number of network membership vectors to compare and N is the
    number of residues. }
  \item{comm.color}{ Either a numeric object or a character object with
    a color name for the plot. }
  \item{important.members}{ A numeric vector corresponding to the
    residue numbers the function should highlight. }
  \item{plot}{ Logical. Whether to make a plot
    representation of the output }
  \item{sse}{ A secondary structure object such as obtained from
    \sQuote{dssp} function to add a secondary structure representation
    to the plot. }
  \item{dash.lines}{ Logical, whether to draw dash lines to separate
    each row in the plot. }
}
\details{
  This function compares the community residue memberships obtained with
  different methods and/or parameters focusing on a specified set of
  residues.
  It shows all the residues that belong to the biggest community
  in the specified residue subset. 
}
\value{
  A numeric binary matrix MxN, where M is the number of community
  clustering methods compared and N the number of nodes in the network
  (i.e. the number of residues in the protein).
  In addition, a plot with all the residues in the biggest community of
  the specified residues can be created.
}
\note{
  If you are providing only a vector membership, please remember to use
  the option drop=FALSE when generating it, otherwise
  \sQuote{membership.comparison} will give an error.
  }
\references{
}
\author{ Guido and Barry }

\seealso{ \code{\link{cna}}, \code{\link{dssp}}}
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

# calculate the secondary structure
# not working for this example
#sse <- dssp(pdb)

# Build, and betweeness cluster, a first network graph 
net1 <- cna(cij)

# Plot the first network
xy1 <- plot.cna(net1)

# Build, and random walk cluster, a second network graph
net2 <- cna(cij,cluster.method="walk")

# Plot the second network
xy2 <- plot.cna(net2)

# Build the membership matrix
memb.matrix <- rbind(net1$raw.communities$membership,
net2$raw.communities$membership)
rownames(memb.matrix) <- c("btwn.nocmap.cij0.4","walk.nocmap.cij0.4")

# Compare the community membership in a specific protein region
membership.comparison(memb.matrix, important.members=c(10:20,100:140), plot=TRUE)

}
\keyword{ utility }