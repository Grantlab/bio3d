

sse.color2 <- function(x, atom.sel=NULL, col.coil="gray", col.helix="purple", col.sheet="yellow") {
	##-- Define a color vector with helix and sheet 
	##    annotations taken from input PDB file.
	##      ToDo: - if no $helix and $sheet defined take 
	## 		        annotation from dssp() or stride()
	##
	## MV inside view.pdb() function
	##
	
	if(!is.null(atom.sel)) {
		x <- trim.pdb(x, atom.sel)
	}
	
	resno <- x$atom$resno
	col <- rep(col.coil, length(resno))

	if(is.null(x$helix) && is.null(x$sheet)) {
		## No SSE defined in PDB calling dssp()
		####x <- dssp(x)
		if(is.null(x$helix) && is.null(x$sheet)) {
			## Still no SSE return coil color
			return(col)
		}
	}

	if(!is.null(x$helix$start)) {
		h.resno <- unbound(x$helix$start, x$helix$end)
		col[(resno %in% h.resno)] = col.helix
	}

	if(!is.null(x$sheet$start)) {
		e.resno <- unbound(x$sheet$start, x$sheet$end)
		col[(resno %in% e.resno)] = col.sheet
	}

	## N.B. The above will not work with multi chain input
	##  i.e. PDB files with overlapping resno, in which 
	##   case we need to match both chain & resno 
	##   - See trim.pdb() for example code
	
	return(col)
}

view <- function(...)
  UseMethod("view")


view.pdb2 <- function(pdb, type="default", atom.sel=NULL, col=NULL, cna=NULL, ...) {
   ##-- Wrapper for visualize() to view larger PDBs the way Barry 
   ##    likes to see them most often.
   ##      To Do - Check validity on "atom.sel" input 
   ##            - Check validity of "col" input and add "keyword" color types
   ##            - Add input more args, e.g. lwd=NULL, lwd.ca=3, lwd.nca=1 etc.
   ##				  - Add cna code.
   ##
   ##     N.B. In general this is still a quick and dirty prototype with 
   ##          no consideration of efficiency and little error checking. 
   ##


   type.options <- c("default", "calpha", "back", "protein", "all")
   type <- match.arg(type, type.options)

   ##- Check on 'col' input 
   ##   This section of code could be improved
   ##   'col' could be a vector of colors or a "keyword" 
   ##    e.g. "sse", "index", "atom", etc. ) 
   if(is.null(col)) {
   	if(type=="all") {
   		## Color by index
   		col <- vec2color(1:nrow(pdb$atom))
   	} else {
   		## Color by secondary structure
   		col <- sse.color2(pdb)
   	}
   } else {
   	if(length(col) == 1) {
   		col=rep(col, nrow(pdb$atom))
   	}
   	if(length(col) != nrow(pdb$atom)) {
   		stop("Length of input color vector, 'col' does not match natom PDB")
   	}
   }

   if(!is.null(atom.sel)) {
   	pdb <- trim.pdb(pdb, atom.sel)
   	col <- col[atom.sel$atom]
   }


   if(type=="default") { 
   	## Calpha trace plus sidechains and ligand
      ca.sel   <- atom.select(pdb, "calpha", verbose = FALSE)
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
      back.sel <- atom.select(pdb, "back", verbose = FALSE)
      lig.sel  <- atom.select(pdb, "ligand", verbose = FALSE)
		side.sel <- combine.select(prot.sel, back.sel, operator="NOT", verbose = FALSE)
		side.sel <- combine.select(side.sel, ca.sel, operator="OR", verbose = FALSE)
  		
  		## Bonds
  		if(!is.null(pdb$con)){
      	connectivity(pdb) <- connectivity(pdb)
      }

      ## Ligand
      if(length(lig.sel$atom) != 0){
        	visualize(trim.pdb(pdb, lig.sel), con=FALSE, type="s")
      }

      ## Sidechain
      if(length(side.sel$atom) != 0) {
      	visualize(trim.pdb(pdb, side.sel), con = FALSE, add = TRUE, type = "l", col = "gray", lwd = 1)
		}

		## Calpha
      if(length(ca.sel$atom) != 0) {
      	ca.pdb <- trim.pdb(pdb, ca.sel)
      	connectivity(ca.pdb) <- calpha.connectivity(ca.pdb)

  			col <- col[ca.sel$atom]

   		visualize(ca.pdb, con=FALSE, col=col, add = TRUE, type = "l", lwd=3, ...)
   	}
    } 


   if(type=="calpha") {
   	## SSE colored C-alpha trace
   	ca.sel <- atom.select(pdb,"calpha", verbose=FALSE)
   	if(length(ca.sel$atom) != 0) {
     		ca.pdb <- trim.pdb(pdb, ca.sel)
     		connectivity(ca.pdb) <- calpha.connectivity(ca.pdb)

  			col <- col[ca.sel$atom]

   		visualize(ca.pdb, con=FALSE, col=col, type = "l", lwd=3, ...)
   	}
   }
   

   if(type=="protein") {
   	## Just protein
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
   	prot.pdb <- trim.pdb(pdb, prot.sel)
	  	connectivity(prot.pdb) <- connectivity(prot.pdb)

		col <- col[prot.sel$atom]

     	if(length(prot.pdb$xyz)!=0){
       	visualize(prot.pdb, con = FALSE, col=col, ...) ## col=col ?
     	}
   }

   
   if(type=="back") {
   	back.sel <- atom.select(pdb,"back", verbose=FALSE)
   	back.pdb <- trim.pdb(pdb, back.sel)
	  	connectivity(back.pdb) <- connectivity(pdb) ## <--- Using back.pdb here gives strange result!!!

  		col <- col[back.sel$atom] 
 
   	if(length(back.pdb$xyz)!=0) {
   		visualize(back.pdb, con=FALSE, col=col, ...)
   	}
   }

   if(type=="all") {
   	## STILL A BUG HERE IN visualize() FOR LIGAND!! - Need to fix
  		### visualize(pdb, con=FALSE, col=col,...) ###

  		## Tmp fix draw ligand separately without lines
  		lig.sel  <- atom.select(pdb, "ligand", verbose = FALSE)
      if(length(lig.sel$atom) != 0){
        	visualize(trim.pdb(pdb, lig.sel), con=FALSE, type="s")
      }

  		prot.pdb <- trim.pdb(pdb, atom.select(pdb,"protein", verbose=FALSE))
	  	connectivity(prot.pdb) <- connectivity(prot.pdb)

		col <- col[prot.sel$atom]

     	if(length(prot.pdb$xyz)!=0){
       	visualize(prot.pdb, con = FALSE, add=TRUE, col=col, ...) ## col=col ?
     	}

  	}
   #if(!is.null(cna)) { 
   #		visualize.cna(cna, pdb, ...) ## need to think more about this
   #}
}

## Need centre=FALSE
view.pdbs <- function(x, type=1, col=NULL, add=FALSE, ...) {
	##-- Wrapper to visualize() for multiple structures

	as.xyz <- function(x, nrow=1, ncol=length(x), byrow=TRUE) {
		y <- matrix(as.numeric(x), nrow=nrow, ncol=ncol, byrow=byrow)
		class(y)="xyz"
		return(y)
	}

	vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
		##-- Define a color scale from a numeric vector
		##     To Do - make independent of classInit package?
		##   (I think Julien posted a new version but I cant find it currently)
		require(classInt)
		return( findColours(classIntervals(vec, n=n, style="equal"), pal) )
	}

	if(class(x) == "3dalign") {
		gap.ind <- is.gap(x$ali)
		x <- x$xyz
	} else {
		###gap.ind <- is.na( x[seq(1, to=length(x), by=3)] ) ##<-- Wrong !!
		gap.ind <- NULL
	}
	nstru <- nrow(x)
	npos  <- (ncol(x)/3)

	## -- The 'type' argument is for trying to sort out 'col' color specification 
	##     for different purposes.  Note. 'col' input could be: 
	##	    1. a) single element vector to be applied to all structures, 
	##         b) a multiple element vector with a color per structure, 
	##      2. a) a vector with a color per atom, or
	##         b) a matrix with a column per atom position and row per structure
	## 
	##     This is specified by 'type=1' or 'type=2' 
	##     Eventually we want to be a bit smarter and remove the need for the 'type' argument 


	## Sort out color options with the aid of type argument
	if(type==1) {
		## Option No. #1 above
		if( is.null(col) ) {
			col <- as.matrix( vmd.colors(nstru) )
		} else {
			if(length(col)==1) {
				col <- as.matrix( rep(col, nstru) )
			}
			if(length(col) != nstru) {
				stop("For type=1: Color vector should be the same length as the number of structures")
			} else {
				col <- as.matrix( col )
			}
		}
	}
	if(type==2) {
		## Option No. #2 above
		if( is.null(col) ) {
			col <- vec2color(1:npos)
			col <- matrix(rep(col, nstru),ncol=npos,nrow=nstru, byrow=TRUE)
		} else {
			if( is.null(nrow(col)) ) {
				## We have an input 'col' vector we want to apply to all structures
				if(length(col) == npos) {
					cat("IN HERE\n\n")
					col <- matrix(rep(col, nstru),ncol=npos,nrow=nstru, byrow=TRUE)
					cat(dim(col))
				}
				else {
					stop("For type=2: Color vector should be same length as ncol pdbs$ali")
					## unclear what the user might want here...
				}
			} else {
				## we have a color matrix
				if(dim(col) != dim(x)) {
					stop("For type=2: Color matrix should be same dim as pdbs$ali")
					## again unclear what the user might want here...
				}
			}
		}
		## Mark gap/missing positions in color matrix for later exclusion
		if( !is.null(gap.ind) )
			col[gap.ind]=NA  ## <---- Trouble here if input is not a pdbs object!!
		##cat( paste(dim(col), collapse="  x  " ), "\n" )
	}

	for(i in 1:nstru) {
		xt = as.xyz( na.omit(x[i,]) )
		xcol = na.omit(col[i,])
		## cat( paste( "  ** length x:", (length(xt)/3), "  length col:", length(xcol),"\n") )
		if(i==1) {
			visualize(xt, con=calpha.connectivity(xt), add=add, col=xcol, ...)
		} else {
			visualize(xt, con=calpha.connectivity(xt), add=TRUE, col=xcol, ...)		
		}
	}
}





stop("dont source more")




##
##- 1. Read and view a PDB file in a useful way 
##
pdb=read.pdb("4q21")

view.pdb2(pdb)
view.pdb2(pdb, "calpha")
view.pdb2(pdb, "back")
view.pdb2(pdb, "all")

view.pdb2(pdb, "trace", col=vec2color(pdb$atom$b))

view.pdb2(pdb, "all", col=vec2color(pdb$atom$b))



##- 2. View a structural fit from a pdbs object!!
data(transducin)
attach(transducin)

view.pdbs(pdbs)   ###<-- lots of 'elesy' related warnings !!! 
view.pdbs(pdbs, col="gray")
view.pdbs(pdbs, col=annotation[,"color"])

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






