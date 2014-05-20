sse.color <- function(x, atom.sel=NULL, col.coil="gray", col.helix="purple", col.sheet="yellow") {
  ##-- Define a color vector with helix and sheet 
  ##    annotations taken from input PDB file.
  ##      ToDo: - if no $helix and $sheet defined take 
  ##            annotation from dssp() or stride()
  ##
  ## MV inside view.pdb() function
  ##
  
#  if(!is.null(atom.sel)) {
#    x <- trim.pdb(x, atom.sel)
#  }
  
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

vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
  col <- colorRampPalette(pal)(n)
  vec.cut <- cut(vec, seq(min(vec), max(vec), length.out=n), include.lowest = TRUE)
  levels(vec.cut) <- 1:length(col)
  col <- col[vec.cut]
  return(col)
}

view <- function(...)
  UseMethod("view")

view.pdb <- function(pdb, type="default", atom.sel=NULL, col=NULL, cna=NULL,
                     elety.custom = atom.index, ...) {
   ##-- Wrapper for visualize() to view larger PDBs the way Barry 
   ##    likes to see them most often.
   ##      To Do - Check validity on "atom.sel" input 
   ##            - Check validity of "col" input and add "keyword" color types
   ##            - Add input more args, e.g. lwd=NULL, lwd.ca=3, lwd.nca=1 etc.
   ##         - Add cna code.
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
      col <- sse.color(pdb)
    }
   } else {
    if(length(col) == 1) {
      col=rep(col, nrow(pdb$atom))
    }
    if(length(col) != nrow(pdb$atom)) {
      stop("Length of input color vector, 'col' does not match natom PDB")
    }
   }

   if(!all(are.symb(pdb$atom$elesy)))
     pdb$atom$elesy <- atom2ele(pdb$atom$elesy, elety.custom)
   
   if(!is.null(atom.sel)) {
    pdb <- trim.pdb(pdb, atom.sel)
    col <- col[atom.sel$atom]
   }

   ## Bonds
   if(is.null(pdb$con)) {
     cat("Computing connectivity from coordinates...\n")
     connectivity(pdb) <- connectivity(pdb)
   }
   
   if(type=="default") { 
      ## Calpha trace plus sidechains and ligand
      ca.sel   <- atom.select(pdb, "calpha", verbose = FALSE)
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
      back.sel <- atom.select(pdb, "back", verbose = FALSE)
      lig.sel  <- atom.select(pdb, "ligand", verbose = FALSE)
      side.sel <- combine.sel(prot.sel, back.sel, op="NOT", verbose = FALSE)
      side.sel <- combine.sel(side.sel, ca.sel, op="OR", verbose = FALSE)
      
      ## Ligand
      if(length(lig.sel$atom) != 0){        
          visualize(trim.pdb(pdb, lig.sel), con=FALSE, type="s", centre=FALSE)
      }

      ## Sidechain
      if(length(side.sel$atom) != 0) {
        visualize(trim.pdb(pdb, side.sel), con = FALSE, add = TRUE,
                  type = "l", col = "gray", lwd = 1, centre=FALSE)        
      }

      ## Calpha
      if(length(ca.sel$atom) != 0) {
        ca.pdb <- trim.pdb(pdb, ca.sel)
        connectivity(ca.pdb) <- calpha.connectivity(ca.pdb)

        col <- col[ca.sel$atom]

        visualize(ca.pdb, con=FALSE, col=col, add = TRUE,
                  type = "l", lwd=3, centre=FALSE, ...)
      }
    } 


   if(type=="calpha") {
    ## SSE colored C-alpha trace
    ca.sel <- atom.select(pdb,"calpha", verbose=FALSE)
    if(length(ca.sel$atom) != 0) {
        ca.pdb <- trim.pdb(pdb, ca.sel)
        connectivity(ca.pdb) <- calpha.connectivity(ca.pdb)

        col <- col[ca.sel$atom]

      visualize(ca.pdb, con=FALSE, col=col, type = "l", lwd=3, centre=FALSE, ...)
    }
   }
   

   if(type=="protein") {
    ## Just protein
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
      prot.pdb <- trim.pdb(pdb, prot.sel)
#       connectivity(prot.pdb) <- connectivity(prot.pdb)

      col <- col[prot.sel$atom]

      if(length(prot.pdb$xyz)!=0){
        visualize(prot.pdb, con = FALSE, col=col, centre=FALSE, ...) ## col=col ?
      }
   }

   
   if(type=="back") {
    back.sel <- atom.select(pdb,"back", verbose=FALSE)
    back.pdb <- trim.pdb(pdb, back.sel)
#      connectivity(back.pdb) <- connectivity(back.pdb) ## <--- Using back.pdb here gives strange result!!!

      col <- col[back.sel$atom] 
 
    if(length(back.pdb$xyz)!=0) {
      visualize(back.pdb, con=FALSE, col=col, centre=FALSE, ...)
    }
   }

   if(type=="all") {
     ## type == all, is it everything but water?
      nowat.sel <- atom.select(pdb, "notwater", verbose = FALSE)
      nowat.pdb <- trim.pdb(pdb, nowat.sel)
      connectivity(nowat.pdb) <- connectivity(nowat.pdb)

      col <- col[nowat.sel$atom]

      if(length(nowat.pdb$xyz)!=0){
        visualize(nowat.pdb, con = FALSE, col=col, centre=FALSE, ...)
      }

    }
   #if(!is.null(cna)) { 
   #    visualize.cna(cna, pdb, ...) ## need to think more about this
   #}
}

view.character <- function(file, type="default", atom.sel=NULL, col=NULL, cna=NULL, ...) {
  x <- read.pdb(file)
  view.pdb(pdb=x, type=type, atom.sel=atom.sel, col=col, cna=cna, ...)
}

## We better add a method for class '3dalign' here, I think.
view.3dalign <- function(x, type=1, col=NULL, add=FALSE, ...) {
  ##-- Wrapper to visualize() for multiple structures
  
  as.xyz <- function(x, ...) {
    x <- matrix(x,nrow=3)
    x <- matrix(x[,!is.na(x[1,])], nrow=1)
    class(x) <- "xyz"
    return(x)
  }
  
  ## 3dalign object should contains an 'xyz' object, not a old style vector!!
  xyz.list <- split(x$xyz, 1:nrow(x$xyz))
  xyz.list <- lapply(xyz.list, as.xyz)
  con.list <- lapply(xyz.list, calpha.connectivity)
  
  ## -- The 'type' argument is for trying to sort out 'col' color specification 
  ##     for different purposes.  Note. 'col' input could be: 
  ##      1. a) single element vector to be applied to all structures, 
  ##         b) a multiple element vector with a color per structure, 
  ##      2. a) a vector with a color per atom, or
  ##         b) a matrix with a column per atom position and row per structure
  ## 
  ##     This is specified by 'type=1' or 'type=2' 
  ##     Eventually we want to be a bit smarter and remove the need for the 'type' argument 
  
  nstru <- length(xyz.list)
  xyz.lengths <- sapply(xyz.list, length)/3
  npos  <- ncol(x$resno)
  ## Sort out color options with the aid of type argument
  if(type==1) {
    ## Option No. #1 above
    if( is.null(col) ) {
      col.list <- vmd.colors(nstru)
    } else {
      if(length(col)==1) {
        col.list <- replicate(col, nstru, simplify = FALSE)
      }
      if(length(col) != nstru) {
        stop("For type=1: Color vector should be the same length as the number of structures")
      } else {
        col.list <- as.list( col )
      }
    }
  }
  if(type==2) {
    ## Option No. #2 above
    if( is.null(col) ) {
      col.list <- lapply(xyz.lengths, function(n) vec2color(1:n))
    } else {
      if(is.list(col)) {
        col.lengths <- lapply(col, length)
        if(any(col.length != xyz.lengths)) {
          stop("When 'col' is a list, each component length must match the number of atoms in each model")
        }
        col.list <- col
      }
      if( is.null(nrow(col)) ) {
        ## We have an input 'col' vector we want to apply to all structures
        if(length(col) == npos) {
#           cat("IN HERE\n\n")
          col    <- replicate(nstru, col, simplify = FALSE)          
          are.na <- lapply(split(x$resno, 1:nstru), is.na)
          col.list <- mapply(function(col, M) return(col[!M]), col, are.na)
#           cat(dim(col))
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
        } else {
          col.list <- split(col, 1:nstru)
          are.na <- lapply(split(x$resno, 1:nstru), is.na)
          col.list <- mapply(function(col, M) return(col[!M]), col.list, are.na)
        }
      }
    }
  }

  elesy <- rep("C", length(xyz.list[[1]])/3)
  visualize(xyz.list[[1]], elesy = elesy, con = con.list[[1]],
            col = col.list[[1]], add = FALSE, centre=FALSE, ...)
  if(length(xyz.list) > 1) {
    useless <- mapply(function(x, con, col){
            elesy <- rep("C", length(x)/3)
            visualize(x, elesy = elesy, con = con, col = col, add = TRUE, centre=FALSE, ...)
           }, xyz.list, con.list, col.list)
  }
}

view.pdbs <- function(x, type=1, col=NULL, add=FALSE, ...) {
  ##-- Wrapper to visualize() for multiple structures

  as.xyz <- function(x, nrow=1, ncol=length(x), byrow=TRUE) {
    y <- matrix(as.numeric(x), nrow=nrow, ncol=ncol, byrow=byrow)
    class(y)="xyz"
    return(y)
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
  ##      1. a) single element vector to be applied to all structures, 
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
      visualize(xt, con=calpha.connectivity(xt), add=add, col=xcol, centre=FALSE, ...)
    } else {
      visualize(xt, con=calpha.connectivity(xt), add=TRUE, col=xcol, centre=FALSE, ...)   
    }
  }
}


