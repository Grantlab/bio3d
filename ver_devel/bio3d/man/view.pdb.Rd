\name{view.pdb}
\alias{view.pdb}
\alias{view.xyz}
\alias{view.pdbs}
\alias{view.character}
\alias{view}
\alias{view.default}
\alias{vec2color}
\title{
  View Molecular Structures 
}
\description{
  Interactive 3D visualization of biomolecular structures from various sources 
  within the Bio3D package. 
}
\usage{
view(...)

\method{view}{pdb}(pdb, as = "all", col = "atom", add = FALSE, elety.custom = NULL, ...)

\method{view}{xyz}(xyz, col = "index", add = FALSE, elesy = NULL, elety.custom = NULL, maxframes = 100, ...)

\method{view}{pdbs}(x, col = "index", add = FALSE, elesy=rep("C", ncol(x$xyz)/3), ...)

vec2color(vec, pal = c("blue", "green", "red"), n = 30)
}
\arguments{
  \item{pdb}{ A PDB object as obtained from function read.pdb.}
  \item{xyz}{ A numeric matrix of xyz coordinates with a frame/structure per row 
    and a Cartesian coordinate per column.}
  \item{x}{ A multi-structure PDBs object as obtained from a number of functions 
    including read.fasta.pdb and pdbaln functions.}
  \item{as}{ A single element character vector that sets the view/render mode 
    determining what and how things are drawn. 
    Current options include \sQuote{all}, \sQuote{overview}, \sQuote{nwat}, 
    \sQuote{protein}, \sQuote{calpha}, \sQuote{back} and \sQuote{ligand} - see 
    examples below for details.}
  \item{col}{ A vector or matrix of colors. The single keywords \sQuote{sse}, 
    \sQuote{index}, \sQuote{atom}, \sQuote{b} and \sQuote{residue} are also valid 
    for \code{\link{view.pdb}}. The keywords \sQuote{index}, \sQuote{atom}, 
    \sQuote{frame}, \sQuote{structure}, \sQuote{gaps} and \sQuote{sse} are valid 
    for \code{\link{view.xyz}} and \code{\link{view.pdbs}} functions. }
  \item{add}{ Logical, if TRUE the requested 3D representation will be added
    to a existing scene. Otherwise a new scene window will be created. }
  \item{elety.custom}{ Used for atom to element type mapping and to
    obtain radii and colors for non standard atoms. }  
  \item{elesy}{ A character vector containing atomic symbols for atoms. }
  \item{maxframes}{ A single element numeric vector specifying the maximum number  
    of frames to be rendered. Note that rendering large trajectories is not the 
    purpose of this function. Rather users should consider VMD and PyMOL for these 
    purposes. } 
  \item{\dots}{ Extra arguments for \code{\link{visualize}}, \code{\link{sse.vector}} 
    and \code{\link{connectivity}} functions. Please see these functions for full 
    argument lists. }
  \item{vec}{ A numeric vector for color mapping. }
  \item{pal}{ A vector of three colors for color pallet specification. }
  \item{n}{ Number of desired color range mappings. }
}
\details{
  The view family of functions generate interactive 3D visualizations of 
  biomolecular structure and dynamics data from various sources within the Bio3D 
  package.  

  The main functions \sQuote{view.pdb} and \sQuote{view.xyz} are wrappers 
  for the \sQuote{visualize.xyz} function. Specialized methods for \sQuote{NMA} 
  and \sQuote{PCA} result visualization are also available along with extensive 
  customization options. 

  Customization options include the ability to set standard graphical parameters 
  such as \sQuote{lwd} for line widths, \sQuote{radii} for atom sphere radii and 
  \sQuote{type} for rendering as points, lines or spheres.

  Extensive \sQuote{col} coloring options are also available. These include coloring 
  by secondary structure (\code{col="sse"}), atom type (\code{col="atom"}), B-factor 
  (\code{col="b"}), residue type (\code{col="residue"}), frame number 
  (\code{col="frame"}), and structure or model number (\code{col="structure"}).

  A matrix or vector of colors matching the input number of atoms and frames can 
  also be used for coloring. Note that for multi-frame input to \sQuote{view.xyz} 
  one can still provide a simple vector of colors but its length should match 
  either the number of atoms (i.e. \code{ncol(xyz)/3}) for atom specific coloring) 
  or number of frames (i.e. \code{nrow(xyz)} for frame specific coloring).

  The \sQuote{as} keyword input argument to \sQuote{view.pdb} currently implements 
  a number of useful representation styles. Current options include: 
  \code{as="all"}, \code{as="overview"}, \code{as="nwat"}, \code{as="protein"}, 
  \code{as="calpha"}, \code{as="back"} and \sQuote{as="ligand"}. See the examples 
  below for details.

  For details of additional settings and options please see the documentation 
  for \code{\link{visualize.xyz}}, \code{\link{connectivity}}, \code{\link{sse.vector}},
  \code{rgl} and \code{par3d}.
}
\value{
  Called for their effect.
}
\references{
  Skjaerven, L. et al. (2014) \emph{BMC Bioinformatics} \bold{15}, 399.
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Barry Grant and Julien Ide }
\seealso{
  \code{\link{visualize}}, \code{\link{pymol}}, 
  \code{rgl.cur}, \code{par3d}, \code{\link{atom.select}},
  \code{\link{vmd.cna}}, \code{\link{vmd.colors}},
  \code{\link{read.pdb}}, \code{\link{write.pdb}},
  \code{\link{trim.pdb}}, \code{\link{cat.pdb}}, 
  \code{\link{read.fasta.pdb}}, \code{\link{as.pdb}},
  \code{\link{read.dcd}}, \code{\link{read.ncdf}}, 
  \code{\link{read.fasta.pdb}}, \code{\link{mktrj.pca}}, 
  \code{\link{mktrj.nma}}
}
\examples{

\dontrun{
##
##- 1. Read and view a PDB file in a useful way 
##
pdb=read.pdb("4q21")

view.pdb(pdb)

## 1. view.pdb examples
view(pdb, "overview")

view(pdb, "calpha", col="sse")
view(pdb, "ligand", add=TRUE)

view(pdb, "back")
view(pdb, "all", col="index")
view(pdb, "calpha", col="b")
view(pdb, "calpha", col="red", lwd=0.5)


##- 2. View a structural fit from a \sQuote{pdbs} class object
data(transducin)
attach(transducin)

view(pdbs)
view(pdbs, col="sse")

view(pdbs, col=annotation[,"color"])
view(pdbs, col="gray")

## define a color scale for RMSF vector
rf=vec2color(rmsf(pdbs$xyz))
view(pdbs, col=rf)



##- 3. View the results of PCA on this structure set (i.e. a \sQuote{xyz} class object)
example(pca.xyz) ## Press RTN.

a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb")
view.xyz(a, col="gray")

# Use d.cut option to increase C-alpha to C-alpha 'trace/bonding' distance if required
view.xyz(a, col=vec2color(rmsf(a)), d.cut=6) 



## Use add=TRUE to add to previous view
view(pdbs)
view(a, col="#808080", add=T)



##- 4. View the results of NMA (i.e. a \sQuote{xyz} class object)
modes <- nma( read.pdb("1hel") )
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")
view.xyz(m7, col=vec2color(rmsf(m7)))  


##- 5. View the results of CNA (i.e. a \sQuote{cna} class object)
#example(plot.cna)
#visualize.cna(net, pdb, xyz.axes=F)
##=> cant turn axis off curently. To be updated as a view function...



##- 6. Simple subregion highlighting
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
pdb <- read.pdb(pdbfile)

## Select and color residues 24 to 27 in both chains
inds <- atom.select(pdb, resno=c(24:27))
mycols <- rep("white", nrow(pdb$atom))
mycols[inds$atom] <- "red"

view(pdb, col=mycols)
view.xyz(pdb$xyz[inds$xyz], col="green", type="s", add=TRUE)


## Motif example highlighting should be easier than the below
##   Ideally allowing updating of the current display with selections

## Lets color motif position
pdb <- read.pdb("5p21")
motif <- "G....GK[ST]"
aa.seq <- pdbseq(pdb)
pos <- motif.find(motif, aa.seq)

aa.seq[pos]

col <- rep("gray", nrow(pdb$atom))
col[ pdb$atom$resno %in% names(aa.seq[pos]) ] = "red"

view(pdb, "calpha", col=col)



##-- Define a color scale for B-factor coloring etc!!
v <- vec2color( pdb$atom$b )
view.pdb(pdb, "overview", col=v)
}
}
\keyword{ utilities }
