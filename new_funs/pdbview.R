
calpha.connectivity <- function(x, d.cut = 4) {
	##-- Quick and dirty Calpha trace connectivity determination
	##    Connects consecutive atoms excluding those more than d.cut 
	##    distance apart
	##     ToDo: - add checks for Calpha only input files
	##           - warn if no connections are found

	if(is.pdb(x)) {
		x <- x$xyz
	}
	natom <- length(x)/3
	con.all <- data.frame(eleno.1 = 1:(natom-1), 
						  eleno.2 = 2:natom)

	a <- matrix(x[atom2xyz(con.all[,1])], ncol=3, byrow=T)
	b <- matrix(x[atom2xyz(con.all[,2])], ncol=3, byrow=T)
	d <- dist.xyz(a, b, all.pairs=FALSE)
	inds <- d < d.cut
	return( con.all[inds,] )
}

sse.color <- function(x, col.coil="gray", col.helix="purple", col.sheet="blue") {
	##-- Define a color vector with helix and sheet 
	##    annotations taken from input PDB file.
	##      ToDo: - if no $helix and $sheet defined take 
	## 		        annotation from dssp() or stride()

	natom <- length(x$xyz)/3
	col <- rep(col.coil, natom)

	h.resno <- unbound(x$helix$start, x$helix$end)
	e.resno <- unbound(x$sheet$start, x$sheet$end)

	h.ind <- which(as.numeric(x$atom[,"resno"]) %in% h.resno)
	e.ind <- which(as.numeric(x$atom[,"resno"]) %in% e.resno)

	col[h.ind]=col.helix
	col[e.ind]=col.sheet
	return(col)
}


vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
	##-- Define a color scale from a numeric vector
	col <- colorRampPalette(pal)(n)
    vec.cut <- cut(vec, seq(min(vec), max(vec), length.out=n),
	               include.lowest = TRUE)
    levels(vec.cut) <- 1:length(col)
    col <- col[vec.cut]
    return(col)
}


pdbview <- function(pdb, sel="default", col=NULL, cna=NULL, ...) {
   ##-- Wrapper for visualize() to view larger PDBs the way Barry 
   ##    likes to see them most often.
   ##      To Do - This is a very quick and dirty prototype with 
   ##              no consideration of efficiency or error checking. 
   ##    For example, we use trim.pdb() rather than just using atom 
   ##    selections and we don't check the 'col' input vector length 
   ##    matches natom or convert it based on 'sel' input).
   ##

   if(sel=="default") { ## this approach is poor!
    ca.pdb <- trim.pdb(pdb, atom.select(pdb,"calpha"))
   	con=calpha.connectivity(ca.pdb)
   	if(is.null(col)) { col=sse.color(ca.pdb, col.coil="#808080") }
   	visualize(ca.pdb, con=con, col=col, ...)

 	## Sidechain
	side.ind <- combine.sel(atom.select(pdb,"protein"), atom.select(pdb,"back"), op="NOT" )
	side.ind <- combine.sel(side.ind, atom.select(pdb,"calpha"), op="OR")
	side.pdb <- trim.pdb(pdb, side.ind)
	visualize(side.pdb, add=TRUE, typ="l", col="gray", lwd=1)

	## Ligand
	lig.pdb <- trim.pdb(pdb, atom.select(pdb,"ligand"))
	visualize(lig.pdb, add=TRUE, typ="s")

   }
   if(sel=="calpha") { ## I like this view!
   	## draw sse colored calpha trace
   	ca.pdb <- trim.pdb(pdb, atom.select(pdb,"calpha"))
   	con=calpha.connectivity(ca.pdb)
   	if(is.null(col)) { col=sse.color(ca.pdb) }
   	visualize(ca.pdb, con=con, col=col, ...)
   } 

   if(sel=="protein") {
   	prot.pdb <- trim.pdb(pdb, atom.select(pdb,"protein"))
   	visualize(prot.pdb, ...)
   }
   if(sel=="back") {
   	back.inds <- combine.sel( atom.select(pdb,"protein"), atom.select(pdb,"back") )
   	prot.pdb <- trim.pdb(pdb, back.inds)
   	if(is.null(col)) { col=sse.color(prot.pdb) }
   	visualize(prot.pdb, col=col, ...)
   }

   #if(!is.null(cna)) { 
   #		visualize.cna(cna, pdb, ...) ## need to think more about this
   #}
}


pdbsview <- function(x, sel, col=NULL, add=FALSE, ...) {
	##-- Wrapper to visualize() for multiple structures
	if(class(x) == "3dalign") {
		x <- pdbs$xyz
	}
	nstru <- nrow(x)

	## -- Abandoned for now trying to sort out color specification in one function
	##    (Note. 'col' input could be: 
	##	    1. single element vector to be applied to all structures, 
	##      2. a vector with a color per structure, 
	##      3. a vector with a color per atom, or
	##      4. a matrix with a column element per position and row per structure
	##     This function can only do #1 & #2. pdbsview2() can do #3 & #4 )

	if( is.null(col) ) {
		col <- vmd.colors(nstru)
	}
	if(length(col)==1) {
		col <- rep(col, nstru)
	}

	for(i in 1:nstru) {
		xt = na.omit(x[i,])
		## if(length(col) == 1)
		if(i==1) {
			visualize(xt, con=calpha.connectivity(xt), add=add, col=col[i], ...)
		###	visualize(xt, con=calpha.connectivity(xt), add=add, col=cols[i,], ...)
		} else {
			visualize(xt, con=calpha.connectivity(xt), add=TRUE, col=col[i], ...)
		###	visualize(xt, con=calpha.connectivity(xt), add=TRUE, col=cols[i,], ...)		
		}
	}
}

pdbsview2 <- function(x, sel, col=NULL, add=FALSE, ...) {
	##-- See pdbsview() for details
	if(class(x) == "3dalign") {
		x <- pdbs$xyz
	}
	nstru <- nrow(x)
	for(i in 1:nstru) {
		xt = na.omit(x[i,])
		## if(length(col) == 1)
		if(i==1) {
			visualize(xt, con=calpha.connectivity(xt), add=add, col=col, ...)
		} else {
			visualize(xt, con=calpha.connectivity(xt), add=TRUE, col=col, ...)
		}
	}
}





##
##- The default output for the visualize function are not that useful for 
##   large structures and also now a little slow with the new data.frame $atom 
##    E.g.
##
pdb=read.pdb("4q21")
visualize(pdb, col=sse.color(pdb)) ## => hard to see fold topology and sse also rotation offset by 'xyz.axes'

## Lets try restricting to Calpha atoms 
ca.pdb <- trim.pdb(pdb, atom.select(pdb,"calpha"))
visualize(ca.pdb, typ="l")               ## => Not useful as no lines drawn
#visualize(ca.pdb, typ="l", safety = 2.7) ## => This gives an Error! ('safety' not passed to connectivity)

ca.con <- connectivity(ca.pdb, safety = 2.7)
visualize(ca.pdb, typ="l", con=ca.con)   ## => Better but some 'bad' connections

## so we need a new connectivity method for Calpha files 
ca.con <- calpha.connectivity(ca.pdb)    ## => This look much better :-)  
visualize(ca.pdb, typ="l", con=ca.con, col=sse.color(ca.pdb))

##--##

##
##-- Overall the above view should be easier to obtain and be able to show C-alpha trace 
##    along with sidechains and ligands. This motivated the following functions as starting 
##    point suggestions for new visualize functions:
##

##- 1. Read and view a PDB file in a useful way 
##
pdb <- read.pdb("1bg2")

pdbview(pdb)            ## Default PDB view 
pdbview(pdb, "calpha")  ## Calpha trace PDB view 
pdbview(pdb, "back")    ## Backbone view
pdbview(pdb, "protein") ## All atom stick view




##- 2. View a structural fit from a pdbs object!!
data(transducin)
attach(transducin)

pdbsview(pdbs, xyz.axes=F, bg.col="black")
##pdbsview2(pdbs, xyz.axes=F, bg.col="black", col=vec2color(rmsf(pdbs$xyz)))



##- 3. View the results of PCA on this structure set
example(pca.xyz)

a <- mktrj.pca(pc.xray, pc=1, file="pc1.pdb")
pdbsview(a, xyz.axes=F, bg.col="black", col="gray")
pdbsview2(a, xyz.axes=F, bg.col="black", col=vec2color(rmsf(a))) ## Cool!!


## can add to previous view
pdbsview(pdbs, xyz.axes=F, bg.col="black")
pdbsview(a, xyz.axes=F, col="#808080", add=T)



##- 4. View the results of NMA
#pdb <- read.pdb("1hel")

modes <- nma(pdb)
m7 <- mktrj.nma(modes, mode=7, file="mode_7.pdb")
pdbsview2(m7, xyz.axes=F, bg.col="black", col=vec2color(rmsf(m7)))
## should pass an argument to calpha.connectivity() to increase 'd.cut' here!


##- 5. View the results of CNA
example(plot.cna)
visualize.cna(net, pdb, xyz.axes=F) ##=> cant turn axis off??? View is offset from center??
## Need to think more about what is most useful here.



##- 6. Simple subregion highlighting should be easier than the below
##      Ideally allowing updating of the current display with selections

## Lets color motif position
motif <- "YRDSKMTRILQDSLGGNCRT"
aa.seq <- pdbseq(pdb)
pos <- motif.find(motif, aa.seq)

natom <- sum(pdb$calpha)
col <- rep("gray", natom)
col[pos]="red"

pdbview(pdb, "calpha", col=col)



##-- Define a color scale for B-factor coloring etc!!

pdb <- read.pdb("1bg2")
v <- as.numeric(pdb$atom[pdb$calpha,"b"])
pdbview(pdb, "calpha", col=vec2color(v))


## lines3d(matrix(ca.pdb$xyz,ncol=3,byrow=TRUE), col=V(net$network)$color)




