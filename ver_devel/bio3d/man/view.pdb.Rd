\name{view.pdb}
\alias{view.pdb}
\alias{view.pdbs}
\alias{view.character}
\alias{vec2color}
\alias{sse.color}
\title{
    View Molecular Structures 
}
\description{
    Interactive 3D visualization of biomolecular structures from various sources within the Bio3D package. 
}
\usage{
view(...)

\method{view}{pdb}(pdb, as = "all", col = "atom", add = FALSE, elety.custom = atom.index, ...)

\method{view}{xyz}(x, type = 1, col = NULL, add = FALSE, ...)

\method{view}{pdbs}(x, type = 1, col = NULL, add = FALSE, ...)

vec2color(vec, pal = c("blue", "green", "red"), n = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdb}{ a PDB object as obtained from read.pdb.}
  \item{as}{ View/render mode that sets what and how things are drawn. 
    Current options include 'all', 'overview', 'nwat', 'protein', 
    'calpha', 'back' and 'ligand' - see examples for details.}
  \item{col}{ A vector or matrix of colors. The single keywords 'sse', 'index',
    'atom', 'b' and 'residue' are also valid for view.pdb. See
    examples. }
  \item{add}{ Logical, if TRUE the requested 3D representation will be added
    to a existing scene. Otherwise a new scene window will be created. }
  \item{elety.custom}{ Used for atom to element type mapping and to
    obtain radii and colors for non standard atoms. }  
  \item{\dots}{ Extra arguments for visualize, sse.vector and connectivity 
    functions. These these functions for full argument list. }
  \item{vec}{ a numeric vector for color mapping. }
  \item{pal}{ color pallet specification. }
  \item{n}{ Number of desired color range mappings. }
  \item{x}{ a numeric coordinate matrix or 3dalign object of pdbs. }
  \item{type}{ A single element numeric vector with values of 1 or 2
  that specifies the Color type rendering for multi-structure viewing 
  - see examples. ToDo: Add option for rmsf coloring and impliment in
  col argument to be consistant with view.pdb also. }
}
\details{
    The view.pdb and view.xyz functions are wrappers for the visualize function.
Basically, they attempt to view larger PDBs in a useful way (i.e. the way we like to see them most often). 
To Do
 - make consistent mode, type and col options across functions. 
 - Write this man page. 
 - Add details of the extensive customization that is possible with visualize.
 - Check validity of "col" input "keywords" for color types
 - Add more input args, e.g. lwd=NULL, lwd.ca=3, lwd.nca=1 etc.
 - Add cna code.
 - Merge view.pdbs and view.xyz codes so only one is needed
 - Add an 'add=FALSE/TRUE' input option so we can add to an existing display.
 
N.B. In general this is still a quick and dirty prototype with 
no consideration of efficiency and little error checking. 

More details to appear here later....
}
\value{
Called for their effect.
}

\examples{
##
##- 1. Read and view a PDB file in a useful way 
##
pdb=read.pdb("4q21")

view.pdb(pdb)

\dontrun{

# 1. view.pdb examples
view(pdb, "overview")

view(pdb, "calpha", col="sse")
view(pdb, "ligand", add=TRUE)

view(pdb, "back")
view(pdb, "all", "index")
view(pdb, "calpha", col="b")
view(pdb, "calpha", "red")


#view(pdb, col="index")
#view(pdb, "calpha", col="index", add=TRUE, lwd=4) 
#view(atom.select(pdb,"ligand",value=TRUE), col="atom", add=TRUE, type="ls")


##- 2. View a structural fit from a '3dalign' class object
data(transducin)
attach(transducin)

view(pdbs)
view(pdbs, type=2)

view(pdbs, col=annotation[,"color"])
view(pdbs, col="gray")

rf=vec2color(rmsf(pdbs$xyz))
view(pdbs, col=rf, type=2)



##- 3. View the results of PCA on this structure set (i.e. a 'xyz' class object)
example(pca.xyz) ## Press RTN.

a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb")
view.xyz(a, col="gray")
view.xyz(a, col=vec2color(rmsf(a)), type=2 ) ## Cool!! 



## can add to previous view
view.xyz(pdbs)
view.xyz(a, col="#808080", add=T)



##- 4. View the results of NMA (i.e. a 'xyz' class object)
#pdb <- read.pdb("1hel")

modes <- nma(pdb)
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")
view.xyz(m7, col=vec2color(rmsf(m7)), type=2)  


##- 5. View the results of CNA (i.e. a 'cna' class object)
#example(plot.cna)
#visualize.cna(net, pdb, xyz.axes=F) ##=> cant turn axis off??? View is offset from center??
## Need to think more about what is most useful here.



##- 6. Simple subregion highlighting should be easier than the below
##      Ideally allowing updating of the current display with selections

## Lets color motif position
motif <- "G....GK[ST]"
aa.seq <- pdbseq(pdb)
pos <- motif.find(motif, aa.seq)

aa.seq[pos]

col <- rep("gray", nrow(pdb$atom))
col[ pdb$atom$resno %in% names(aa.seq[pos]) ] = "red"

view(pdb, "calpha", col=col)



##-- Define a color scale for B-factor coloring etc!!
v <- vec2color( pdb$atom$b )
view.pdb(pdb, "trace", col=v)
}
}
\keyword{ utilities }
