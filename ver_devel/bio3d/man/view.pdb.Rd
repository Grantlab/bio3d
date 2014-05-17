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
view.pdb(pdb, type = "default", atom.sel = NULL, col = NULL, cna = NULL, ...)
view.pdbs(x, type = 1, col = NULL, add = FALSE, ...)
vec2color(vec, pal = c("blue", "green", "red"), n = 30)
sse.color(x, atom.sel=NULL, col.coil="gray", col.helix="purple", col.sheet="yellow") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdb}{ a PDB object as obtained from read.pdb.}
  \item{type}{ Render type - see examples.}
  \item{atom.sel}{ A selection object from atom.select}
  \item{col}{ A vector or matrix of colors.}
  \item{cna}{ A protein structure network object. Not yet implemented. }
  \item{\dots}{ Extra arguments for visualize functions. }
  \item{vec}{ a numeric vector for color mapping. }
  \item{pal}{ color pallet specification. }
  \item{n}{ Number of desired color range mappings. }
  \item{x}{ a numeric coordinate matix or 3dalign object of pdbs. }
}
\details{
    The view.pdb and view.pdbs functions are wrappers for the visualize function.
Basically, they attempt to view larger PDBs the way Barry likes to see them most often.
To Do
 - Write this man page
 - Check validity on "atom.sel" input 
 - Check validity of "col" input and add "keyword" color types
 - Add input more args, e.g. lwd=NULL, lwd.ca=3, lwd.nca=1 etc.
 - Add cna code.

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

view.pdb(pdb, "calpha")
view.pdb(pdb, "back")
view.pdb(pdb, "all")

#view.pdb(pdb, "trace", col=vec2color(pdb$atom$b))

view.pdb(pdb, "all", col=vec2color(pdb$atom$b))



##- 2. View a structural fit from a pdbs object!!
data(transducin)
attach(transducin)

# view.pdbs(pdbs)   ###<-- lots of elesy related warnings !!! 
view.pdbs(pdbs, col=annotation[,"color"])
view.pdbs(pdbs)   ###<-- lots of 'elesy' related warnings !!! 
view.pdbs(pdbs, col="gray")

rf=vec2color(rmsf(pdbs$xyz))
view.pdbs(pdbs, col=rf, type=2)



##- 3. View the results of PCA on this structure set
example(pca.xyz) ## Press RTN.

a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb")
view.pdbs(a, col="gray")
view.pdbs(a, col=vec2color(rmsf(a)), type=2 ) ## Cool!! 



## can add to previous view
view.pdbs(pdbs)
view.pdbs(a, col="#808080", add=T)



##- 4. View the results of NMA
#pdb <- read.pdb("1hel")

modes <- nma(pdb)
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")
view.pdbs(m7, col=vec2color(rmsf(m7)), type=2)  


##- 5. View the results of CNA
#example(plot.cna)
#visualize.cna(net, pdb, xyz.axes=F) ##=> cant turn axis off??? View is offset from center??
## Need to think more about what is most useful here.



##- 6. Simple subregion highlighting should be easier than the below
##      Ideally allowing updating of the current display with selections

## Lets color motif position
motif <- "G....GK[ST]"
aa.seq <- pdbseq(pdb)
pos <- motif.find(motif, aa.seq)

natom <- sum(pdb$calpha)
col <- rep("gray", natom)
col[pos]="red"

pdbview(pdb, "trace", col=col)



##-- Define a color scale for B-factor coloring etc!!
v <- vec2color( pdb$atom$b )
view.pdb(pdb, "trace", col=v)
}
}
\keyword{ utilities }
