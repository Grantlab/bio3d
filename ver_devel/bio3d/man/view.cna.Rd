\name{view.cna}
\alias{view.cna}
\title{
  View cna commumnity output in VMD
}
\description{
  This function generates a file that can be superimposed to the protein
  structure in VMD
}
\usage{
view.cna(x, pdb, layout = layout.pdb(pdb, x), 
		col.sphere=NULL, col.lines = "silver", weights = NULL, 
		radius = table(x$communities$membership)/5, alpha = 1, 
		vmdfile = "network.vmd", pdbfile = "network.pdb", 
		launch = FALSE)
}
\arguments{
  \item{x}{A 'cna' class object such as obtained from \sQuote{cna} function/
}
  \item{pdb}{A 'pdb' class object such as obtained from \sQuote{read.pdb} function.
}
  \item{layout}{ A numeric matrix of Nx3, where N is the number of spheres
(i.e. communities) to be drawn.
}
  \item{col.sphere}{ A numeric vector containing the sphere colors.
}
  \item{col.lines}{ A character object specifying the color of the
    edges (default 'silver'). Must use VMD colors names.
}
  \item{weights}{ A numeric vector specifying the edge width. Default is
    taken from E(x$clustered.network)$weight.
}
  \item{radius}{ A numeric vector containing the sphere radii. Default
    is taken from the number of community members divided by 5.
}
  \item{alpha}{ A numeric element specifying VMD alpha
    parameter. Default is set to 1.
}
  \item{vmdfile}{ A character element specifying the drawing file name to be
    loaded in VMD. The residue chain is assigned according to the
    community number (it makes easier to color the communities).
}
  \item{pdbfile}{ A character element specifying the pdb file name to be
    loaded in VMD.
}
  \item{launch}{ Logical. If TRUE, a VMD session will be started with
    the output of \sQuote{view.cna}
}
}
\details{ This function generates a ball and stick representation of the
  community clustering that overlaps the protein structure. It can be
  visualized in VMD. The sphere radii and the edge thickness are
  proportional to the community members and the edge weights, respectively.
}
\value{
Two files are generated as output. A pdb file with the residue chains
assigned according to the community and a text file containing the
drawing commands for the community representation.
}
\references{ Humphrey, W., Dalke, A. and Schulten, K., ``VMD - Visual Molecular Dynamics'' J. Molec. Graphics 1996, 14.1, 33-38.
}
\author{
Barry Grant}

\examples{
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

# View cna
view.cna(net, pdb, layout=layout.pdb(pdb,net), launch=FALSE)
}

\keyword{ utility }
