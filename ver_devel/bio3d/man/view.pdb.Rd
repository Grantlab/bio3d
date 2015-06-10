\name{view.pdb}
\alias{view.pdb}
\alias{view.xyz}
\alias{view.pdbs}
\alias{view.character}
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

\method{view}{pdb}(pdb, as = "all", col = "atom", add = FALSE, elety.custom = atom.index, ...)

\method{view}{xyz}(xyz, col = "index", add = FALSE, elesy = NULL, elety.custom = atom.index, ...)

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
    Current options include \sQuote{all}, \sQuote{overview}, \sQuote{nwat}, \sQuote{protein}, 
    \sQuote{calpha}, \sQuote{back} and \sQuote{ligand} - see examples below for details.}
  \item{col}{ A vector or matrix of colors. The single keywords \sQuote{sse}, \sQuote{index},
    \sQuote{atom}, \sQuote{b} and \sQuote{residue} are also valid for \code{\link{view.pdb}}. The keywords \sQuote{index}, 
    \sQuote{atom}, \sQuote{frame}, \sQuote{structure}, \sQuote{gaps} and \sQuote{sse} are valid for \code{\link{view.xyz}} and 
    \code{\link{view.pdbs}} functions. }
  \item{add}{ Logical, if TRUE the requested 3D representation will be added
    to a existing scene. Otherwise a new scene window will be created. }
  \item{elety.custom}{ Used for atom to element type mapping and to
    obtain radii and colors for non standard atoms. }  
  \item{\dots}{ Extra arguments for \code{\link{visualize}}, \code{\link{sse.vector}} and \code{\link{connectivity}} 
    functions. Please see these functions for full argument lists. }
  \item{vec}{ A numeric vector for color mapping. }
  \item{pal}{ A vector of three colors for color pallet specification. }
  \item{n}{ Number of desired color range mappings. }
}
\details{
  The view() family of functions generate interactive 3D visualizations of 
  biomolecular structure and dynamics data from various sources within the Bio3D 
  package.  The main functions \sQuote{view.pdb} and \sQuote{view.xyz} are wrappers for the 
  \sQuote{visualize.xyz} function. Extensive customization are possible, including 
  standard graphical parameters such as \sQuote{lwd} and more advanced representation 
  styles and colors with the \sQuote{as} and \sQuote{col} keyword options. Full details of these customizations will be added here as these functions mature. For now please see 
  below examples for further details as well as the documentation for \sQuote{visualize.xyz'.

  Current \sQuote{as} options include: \sQuote{all}, \sQuote{overview}, \sQuote{nwat}, \sQuote{protein}, \sQuote{calpha}, 
  \sQuote{back} and \sQuote{ligand'.

  Current \sQuote{col} keyword options for \sQuote{view.pdb} include: \sQuote{sse}, \sQuote{index}, \sQuote{atom}, 
  \sQuote{b} and \sQuote{residue'.  The keywords \sQuote{index}, \sQuote{atom}, \sQuote{frame}, \sQuote{structure}, \sQuote{gaps} and \sQuote{sse} are valid for view.pdbs and view.xyz functions. 
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


##- 2. View a structural fit from a \sQuote{pdbs} class object
data(transducin)
attach(transducin)

view(pdbs)
view(pdbs, col="sse")

view(pdbs, col=annotation[,"color"])
view(pdbs, col="gray")

rf=vec2color(rmsf(pdbs$xyz))
view(pdbs, col=rf)



##- 3. View the results of PCA on this structure set (i.e. a \sQuote{xyz} class object)
example(pca.xyz) ## Press RTN.

a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb")
view.xyz(a, col="gray")
view.xyz(a, col=vec2color(rmsf(a)), d.cut=6) ## Cool!! 



## can add to previous view
view(pdbs)
view(a, col="#808080", add=T)



##- 4. View the results of NMA (i.e. a \sQuote{xyz} class object)
#pdb <- read.pdb("1hel")

modes <- nma(pdb)
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")
view.xyz(m7, col=vec2color(rmsf(m7)))  


##- 5. View the results of CNA (i.e. a \sQuote{cna} class object)
#example(plot.cna)
#visualize.cna(net, pdb, xyz.axes=F) ##=> cant turn axis off??? View is offset from center??
## Need to think more about what is most useful here.



##- 6. Simple subregion highlighting should be easier than the below
##      Ideally allowing updating of the current display with selections

## Lets color motif position
pdb<-read.pdb("5p21")
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
