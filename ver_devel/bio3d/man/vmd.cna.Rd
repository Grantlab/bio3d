\name{vmd}
\alias{vmd}
\alias{vmd.cna}
\alias{vmd.cnapath}
\title{ View CNA Protein Structure Network Community Output in VMD }
\description{
  This function generates a VMD scene file and a PDB file that can be 
  read and rendered by the VMD molecular viewer. Chose \sQuote{color by chain} 
  to see corresponding regions of structure colored by community along with 
  the community protein structure network. 
}
\usage{
vmd(\dots)

\method{vmd}{cna}(x, pdb, layout = layout.cna(x, pdb, k=3), 
		col.sphere=NULL, col.lines = "silver", weights = NULL, 
		radius = table(x$communities$membership)/5, alpha = 1,
		vmdfile = "network.vmd", pdbfile = "network.pdb", 
                full = FALSE, launch = FALSE, exefile=NULL, \dots)

\method{vmd}{cnapath}(x, pdb, out.prefix = "vmd.cnapath", spline = FALSE, 
   colors = c("blue", "red"), launch = FALSE, exefile=NULL, \dots)

}
\arguments{
  \item{x}{A 'cna' or 'cnapath' class object as obtained from functions
    \code{cna} or \code{cnapath}.
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
  \item{full}{ Logical, if TRUE the full all-atom network rather than the
    clustered community network will be drawn. Intra community edges are
    colored according to the community membership, while inter community
    edges are thicker and colored black.
  }
  \item{launch}{ Logical. If TRUE, a VMD session will be started with
    the output of \sQuote{vmd.cna}.
  }
  \item{out.prefix}{ Prefix for the names of output files, \sQuote{vmd.cnapath.vmd}
    and \sQuote{vmd.cnapath.pdb}.
  }
  \item{spline}{ Logical, if TRUE all paths are displayed as spline
    curves. 
  }
  \item{colors}{
    Character vector or integer scalar, define path colors. If a
    character vector, passed to \code{\link{colorRamp}} function to generate the 
    color scales. If an integer, color all paths the same way with VMD color 
    ID equal to the integer. }
  \item{exefile}{ file path to the \sQuote{VMD} program on your system (i.e.
    how is \sQuote{VMD} invoked). If \code{NULL}, use OS-dependent
    default path to the program. }
  \item{\dots}{ additional arguments passed to the function
    \code{\link{colorRamp}} (in \code{vmd.cnapath}).
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

if (!requireNamespace("igraph", quietly = TRUE)) {
   message('Need igraph installed to run this example')
} else {

# Load the correlation network from MD data
attach(hivp)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

# View cna 
vmd.cna(net, pdb, launch=FALSE)
## within VMD set 'coloring method' to 'Chain' and 'Drawing method' to Tube


##-- From NMA
pdb.gdi = read.pdb("1KJY")
pdb.gdi = trim.pdb(pdb.gdi, inds=atom.select(pdb.gdi, chain="A", elety="CA")) 
modes.gdi = nma(pdb.gdi)
cij.gdi = dccm(modes.gdi)
net.gdi = cna(cij.gdi, cutoff.cij=0.35)
#vmd.cna(net.gdi, pdb.gdi, alpha = 0.7, launch=TRUE)

detach(hivp)

}
}
}

\keyword{ utility }
