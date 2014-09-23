\name{view.cna}
\alias{view.cna}
\title{
  View CNA Protein Structure Network Community Output in VMD
}
\description{
  This function generates a VMD scene file and a PDB file that can be 
  read and rendered by the VMD molecular viewer. Chose \sQuote{color by chain} 
  to see corresponding regions of structure colored by community along with 
  the community protein structure network. 
}
\usage{
view.cna(x, pdb, layout = layout.cna(x, pdb, k=3), 
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
  \item{layout}{ A numeric matrix of Nx3 XYZ coordinate matrix, where N is the 
  number of community spheres to be drawn.
  }
  \item{col.sphere}{ A numeric vector containing the sphere colors.
  }
  \item{col.lines}{ A character object specifying the color of the
    edges (default 'silver'). Must use VMD colors names.
  }
  \item{weights}{ A numeric vector specifying the edge width. Default is
    taken from E(x$community.network)$weight.
  }
  \item{radius}{ A numeric vector containing the sphere radii. Default
    is taken from the number of community members divided by 5.
  }
  \item{alpha}{ A single element numeric vector specifying the VMD alpha
    transparency parameter. Default is set to 1.
  }
  \item{vmdfile}{ A character element specifying the output VMD scene file 
    name that will be loaded in VMD.
  }
  \item{pdbfile}{ A character element specifying the output pdb file name to be
    loaded in VMD.
  }
  \item{launch}{ Logical. If TRUE, a VMD session will be started with
    the output of \sQuote{view.cna}.
}
}
\details{ 
  This function generates a scaled sphere (communities) and stick (edges) 
  representation of the community network along with the corresponding protein 
  structure divided into chains, one chain for each community. The sphere radii 
  are proportional to the number of community members and the edge widths correspond 
  to network edge weights.
}
\value{
  Two files are generated as output. A pdb file with the residue chains
  assigned according to the community and a text file containing The
  drawing commands for the community representation.
}
\references{ Humphrey, W., Dalke, A. and Schulten, K., ``VMD - Visual Molecular Dynamics'' J. Molec. Graphics 1996, 14.1, 33-38.
}
\author{
Barry Grant}

\examples{
\dontrun{
  ##-- From MD
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# select residues 24 to 27 and 85 to 90 in both chains
inds <- atom.select(pdb, resno=c(24:27, 85:90), elety="CA")

# lsq fit of trj on pdb
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=inds$xyz, mobile.inds=inds$xyz)

# calculate dynamical cross-correlation matrix
cij <- dccm(xyz)
#plot(cij)

# Build, and betweeness cluster, a network graph
net <- cna(cij)

# View cna 
view.cna(net, pdb, launch=FALSE)
## within VMD set 'coloring method' to 'Chain' and 'Drawing method' to Tube


##-- From NMA
pdb.gdi = read.pdb("1KJY")
pdb.gdi = trim.pdb(pdb.gdi, inds=atom.select(pdb.gdi, chain="A", elety="CA")) 
modes.gdi = nma(pdb.gdi)
cij.gdi = dccm(modes.gdi)
net.gdi = cna(cij.gdi, cutoff.cij=0.35)
view.cna(net.gdi, pdb.gdi, alpha = 0.7, launch=TRUE)
}
}

\keyword{ utility }
