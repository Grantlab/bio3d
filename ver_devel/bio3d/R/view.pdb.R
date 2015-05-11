
vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
  col <- colorRampPalette(pal)(n)
  vec.cut <- cut(vec, seq(min(vec), max(vec), length.out=n), include.lowest = TRUE)
  levels(vec.cut) <- 1:length(col)
  col <- col[vec.cut]
  return(col)
}

view <- function(...)
  UseMethod("view")

view.pdb <- function(pdb, mode="default", atom.sel=NULL, col=NULL, cna=NULL,
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


   mode.options <- c("default", "calpha", "back", "protein", "all")
   mode <- match.arg(mode, mode.options)

   ##- Check on 'col' input
   ##   This section of code could be improved
   ##   'col' could be a vector of colors or a "keyword"
   ##    e.g. "sse", "index", "atom", etc. )
   if(is.null(col)) {
    if(mode=="all") {
      ## Color by index
      col <- vec2color(1:nrow(pdb$atom))
    } else {
      ## Color by secondary structure
      ##col <- sse.color(pdb)
      col <- sse.vector(pdb, coil="gray", helix="purple", sheet="yellow", calpha=FALSE)
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

   if(mode=="default") {
      ## Calpha trace plus sidechains and ligand
      ca.sel   <- atom.select(pdb, "calpha", verbose = FALSE)
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
      back.sel <- atom.select(pdb, "back", verbose = FALSE)
      lig.sel  <- atom.select(pdb, "ligand", verbose = FALSE)
      side.sel <- combine.select(prot.sel, back.sel, operator="NOT", verbose = FALSE)
      side.sel <- combine.select(side.sel, ca.sel, operator="OR", verbose = FALSE)

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


   if(mode=="calpha") {
    ## SSE colored C-alpha trace
    ca.sel <- atom.select(pdb,"calpha", verbose=FALSE)
    if(length(ca.sel$atom) != 0) {
        ca.pdb <- trim.pdb(pdb, ca.sel)
        connectivity(ca.pdb) <- calpha.connectivity(ca.pdb)

        col <- col[ca.sel$atom]

      visualize(ca.pdb, con=FALSE, col=col, type = "l", lwd=3, centre=FALSE, ...)
    }
   }


   if(mode=="protein") {
    ## Just protein
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
      prot.pdb <- trim.pdb(pdb, prot.sel)
#       connectivity(prot.pdb) <- connectivity(prot.pdb)

      col <- col[prot.sel$atom]

      if(length(prot.pdb$xyz)!=0){
        visualize(prot.pdb, con = FALSE, col=col, centre=FALSE, ...) ## col=col ?
      }
   }


   if(mode=="back") {
    back.sel <- atom.select(pdb,"back", verbose=FALSE)
    back.pdb <- trim.pdb(pdb, back.sel)
#      connectivity(back.pdb) <- connectivity(back.pdb) ## <--- Using back.pdb here gives strange result!!!

      col <- col[back.sel$atom]

    if(length(back.pdb$xyz)!=0) {
      visualize(back.pdb, con=FALSE, col=col, centre=FALSE, ...)
    }
   }

   if(mode=="all") {
     ## mode == all, is it everything but water?
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

view.character <- function(file, mode="default", atom.sel=NULL, col=NULL, cna=NULL, ...) {
  ## This function is probably not needed
  ## Useful to visualize quickly a pdb file using a filename
  x <- read.pdb(file)
  view.pdb(pdb=x, mode=mode, atom.sel=atom.sel, col=col, cna=cna, ...)
}

view.xyz <- function(x, type=1, col=NULL, add=FALSE, ...) {
  ##-- Wrapper to visualize() for multiple structures
  ##
  ## ToDo. - Incorporate best bits of view.3dalign()
  ##       - Ask user if more than 300 frames are to be drawn
  ##       - Generally speed up - no need to re-calculate connectivity
  ##          for xyz input (but required for 3dalign input)
  if(!inherits(x, "xyz"))
    stop("'x' must be an object of class 'xyz'")

  as.xyz <- function(x, ...) {
    #     x <- matrix(x,nrow=3)
    x <- matrix(x[,!is.na(x[1,])], nrow=1)
    class(x) <- "xyz"
    return(x)
  }

  nstru <- nrow(x)
  npos  <- ncol(x)/3

  if(nstru > 300) {
    warning( paste("Input 'x' has",nstru, "frames. Only drawing first 300") )
    x <- x[1:300,]
    nstru <- 300
  }

  xyz.list <- split(x, 1:nrow(x))
  xyz.list <- lapply(xyz.list, matrix, nrow = 3)
  are.na.list <- lapply(xyz.list, function(x) return(is.na(x[1,])))
  xyz.list <- lapply(xyz.list, as.xyz)
  xyz.lengths <- sapply(xyz.list, length)/3

  ## Compute the connectivity only once if all the xyz.models have the same number of atoms (no NA)
  model.with.na <- sapply(are.na.list, any)
  if(any(model.with.na)) {
    con.list <- lapply(xyz.list, calpha.connectivity)
  } else {
    con.list <- replicate(nstru, calpha.connectivity(xyz.list[[1]]), simplify = FALSE)
  }

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
          stop("When 'col' is a list, each component length must match the number of atoms/residue in each model")
        }
        col.list <- col
      }
      if( is.null(nrow(col)) ) {
        ## We have an input 'col' vector we want to apply to all structures
        if(length(col) == npos) {
          #           cat("IN HERE\n\n")
          col    <- replicate(nstru, col, simplify = FALSE)
          col.list <- mapply(function(col, M) return(col[!M]), col, are.na.list)
          #           cat(dim(col))
        }
        else {
          stop("For type=2: Color vector should be same length as ncol of your 'xyz' object")
          ## unclear what the user might want here...
        }
      } else {
        ## we have a color matrix
        if(dim(col) != dim(x)) {
          stop("For type=2: Color matrix should be same dim as your 'xyz' object")
          ## again unclear what the user might want here...
        } else {
          col.list <- split(col, 1:nstru)
          col.list <- mapply(function(col, M) return(col[!M]), col.list, are.na.list, SIMPLIFY=FALSE)
        }
      }
    }
  }

  elesy <- rep("C", length(xyz.list[[1]])/3)
  visualize(xyz.list[[1]], elesy = elesy, con = con.list[[1]],
            col = col.list[[1]], add = add, centre=FALSE, ...)
  if(length(xyz.list) > 1) {
    useless <- mapply(function(x, con, col){
      elesy <- rep("C", length(x)/3)
      visualize(x, elesy = elesy, con = con, col = col, add = TRUE, centre=FALSE, ...)
    }, xyz.list, con.list, col.list)
  }
}

view.pdbs <- function(x, type=1, col=NULL, add=FALSE, ...) {
  ##-- Wrapper to visualize() for multiple structures
  ##
  ## ToDo.   Combine/merge with view.xyz() below and then simply
  ##          call view.xyz() within view.pdbs() view.nma()
  ##          view.pca(), view.cna() etc.
  if(!inherits(x, "pdbs"))
    stop("'x', must be an object of class 'pdbs' as obtained from read.fasta.pdb()")
  view.xyz(x$xyz, type = type, col = col, add = add, ...)
}
