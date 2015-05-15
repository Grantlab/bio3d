
vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
  col <- colorRampPalette(pal)(n)
  vec.cut <- cut(vec, seq(min(vec), max(vec), length.out=n), include.lowest = TRUE)
  levels(vec.cut) <- 1:length(col)
  return(col[vec.cut])
}


### Testing snippets
#pdb <- read.pdb("5p21")
#y <- trim.pdb(pdb, atom.select(pdb, "sideCA", verbose = FALSE))
#
#view(pdb)                             ## DULL default! 
#view(pdb, "overview", "sse", add=T)   ## can add on top
#
#view(pdb, col="sse", helix="green", lwd=2, xyz.axes=T) ## can pass args intelligently   
#
#view(pdb, "calpha")                   ## I miss the sse coloring
#view(pdb, "calpha", "sse")            ##  there it is but I have to type more
#
#view(pdb, "calpha", "sse", type="s", add=T)
#view(y, add=TRUE, col="white", lwd=3)
#view(atom.select(pdb,"ligand",value=T), col="atom", add=T, type="ls")
#
#view(pdb, col="index")
#view(pdb, "calpha", col="index", add=T, lwd=4) 
#view(atom.select(pdb,"ligand",value=T), col="atom", add=T, type="ls")

view <- function(...)
  UseMethod("view")


view.pdb <- function(pdb, as="all", col="atom", add=FALSE, 
                     elety.custom = atom.index, ...) {

  ##-- Wrapper for visualize() to view larger PDBs the way Barry
  ##    likes to see them most often.
  ##    Note. We can now use all visualize input args like: lwd, bg.col, etc.
  ##          We can also add to existing view by setting add=T
  ##      To Do           
  ##            - Check for duplicated input args, and take sse cols and 
  ##               also calpha.connectivity d.cut from \dots 
  ##            - Add color by residue type code. 
  ##
  ##     N.B. In general this is still a quick and dirty prototype with
  ##          no consideration of efficiency and little error checking.
  ##




  arg.filter <- function(new.args, FUN=NULL, dots=list(...)) {
    ##-- Simple list filtering for duplicate 
    ##    function input argument removal and validation.
    ##
    ## new.args = The new default args that can be overwritten 
    ##              by those in 'dots' (i.e. user supplied "...") 
    ##               E.G. "new.args=list(col=mydefualtcol, lwd=3)"
    ##
    ## FUN      = Function name from which allowed arguments are checked
    ##              and used to limit output of this function.
    ##              This is typically only required if there are multiple 
    ##              (sub)functions to be called each with other specific  
    ##              things in /dots.
    ##               E.G. allowed=names(formals( mysubfunction2call ))
    ##
    ## dots     = Full user supplied updated values typically 
    ##               this is the extra /dots values, i.e. list(...) 
    ##    
    ##   sse.default <- list(coil="gray", helix="purple", sheet="yellow")
    ##   sse.args <- arg.filter( sse.default, FUN=sse.vector )
    ##   col <- do.call('sse.vector', c(list(pdb=pdb), sse.args) )
    ##
    ## Returns entries of 'dots' updated with those in 'new.args'
    ##   that intersect with allowed FUN input args.

    ans <- c(dots, new.args[!names(new.args) %in% names(dots)])
    if(!is.null(FUN)) { ans <- ans[names(ans) %in% names(formals(FUN))] }
    return(ans)
  }


  ##- Check for possible input 'd.cut' to be used in calpha.connectivity()
  d <- list(...)$d.cut
  if(is.null(d)) { d <- 4 }
  if(d < 4) { stop("Input 'd.cut' should be 4 or more to avoid rendering errors")}


  ##-- Input 'as=' option will set what is actually drawn
  as.options <- c("all", "overview", "nwat", "protein", "calpha", "back", "ligand")
  as <- match.arg(as, as.options)

  ##-- Input 'col=' option will set how drawn things are colored 
  display <- NULL
  if(length(col) == 1) {
    ##- Requested 'col' may be a 'display color keyword' rather than a simple color
    display.options <- c("sse", "index", "atom", "b", "residue")
    display.check <- pmatch(col, display.options, nomatch = NA)
    if(!is.na( display.check )) {
      display <- display.options[display.check]
    }
  } 


  ##- Using 'Display color keywords'
  if( !is.null(display) ){
    ## Do we have any protein and therefore need to consider SSE coloring?
    prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
    input.prot <- (length(prot.sel$atom) > 3)

    ## Color by secondary structure if protein present
    if(display=="sse") {
      if(input.prot) {
          sse.default <- list(coil="gray", helix="purple", sheet="yellow", calpha=FALSE)
          sse.args <- arg.filter( sse.default, FUN=sse.vector )
          col <- do.call('sse.vector', c(list(pdb=pdb), sse.args) )
      } else {
        col <- "gray"
      }
    }
  
    ## Color by index position 
    if(display=="index") { col <- vec2color(1:nrow(pdb$atom)) }

    ## Color by atom type 
    if(display=="atom") { col <- NULL }

    ## Color by numeric vector in b-factor field
    if(display=="b") { col <- vec2color(as.numeric(pdb$atom$b)) }

    ## Color by residue type  
    if(display=="residue") {
      cat("Residue coloring not yet implemented, Sorry!\n")
      col <- NULL
    }

  } else {
    ##- Check on input color option length vs Natom
    if(length(col) != nrow(pdb$atom)) {
      if(!is.null(col) && (length(col) > 1))
        warning("Length of input color vector, 'col', does not match natom PDB: recycling colors")

      col=rep(col, nrow(pdb$atom))
    }
  }

  ##- Atom/Element type check
  if(!all(are.symb(pdb$atom$elesy)))
    pdb$atom$elesy <- atom2ele(pdb$atom$elesy, elety.custom)


  ## Bonds/Connectivity check
  #if(is.null(pdb$con)) {
  #  cat("Computing connectivity from coordinates...\n")
  #  connectivity(pdb) <- connectivity(pdb)
  #}

  ##- Store often used visualize arguments
  vis.default.args <- arg.filter( list(con=FALSE, col=col, add=add, centre=FALSE), FUN=visualize.pdb )


  if(as=="all") {
    ## Here we just use visualize() with our new default args, plus an
    ##  extra catch for calpha-only input (for alternate connectivity)
    if((length(atom.select(pdb,"protein")$atom) == sum(pdb$calpha)) && (sum(pdb$calpha) > 1))
      connectivity(pdb) <- calpha.connectivity(pdb, d.cut=d)

    do.call('visualize.pdb', c(list(pdb=pdb), vis.default.args) )

    ## Plus a light non-protein atom display.
    np <- atom.select(pdb,"notprotein")
    if(length(np$atom) > 0)
      visualize(trim.pdb(pdb, np), type="p", con = FALSE, centre=FALSE, add=TRUE)
  }

  if(as=="overview") {
    ## Calpha trace plus sidechains and ligand
    ca.sel   <- atom.select(pdb, "calpha", verbose = FALSE)
    prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
    back.sel <- atom.select(pdb, "back", verbose = FALSE)
    lig.sel  <- atom.select(pdb, "ligand", verbose = FALSE)
    side.sel  <- atom.select(pdb, "side", verbose = FALSE)

    ## Ligand - hardwired 'type' etc. for now
    if(length(lig.sel$atom) != 0){
      visualize(trim.pdb(pdb, lig.sel), con=FALSE, type="s", centre=FALSE, add=add)
      add <- TRUE
    }

    ## Sidechain - note hardwired 'lwd' and 'type'
    if(length(side.sel$atom) != 0) {
      side.sel <- combine.select(side.sel, ca.sel, operator="OR", verbose = FALSE)
      visualize(trim.pdb(pdb, side.sel), con = FALSE, add = add,
                type = "l", col=col[side.sel$atom], lwd = 1, centre=FALSE)
      add <- TRUE
    }

    ## Calpha Trace
    if(length(ca.sel$atom) != 0) {
      ca.pdb <- trim.pdb(pdb, ca.sel)
      connectivity(ca.pdb) <- calpha.connectivity(ca.pdb, d.cut=d)
      col <- col[ca.sel$atom]

      vis.default <- list(type="l", lwd=4, con=FALSE, col=col, add=add, centre=FALSE) 
      vis.args <- arg.filter( vis.default, FUN=visualize.pdb )
      do.call('visualize.pdb', c(list(pdb=ca.pdb), vis.args) )
     }
  }

  if(as=="calpha") {
    ## C-alpha trace only
    ca.sel <- atom.select(pdb,"calpha", verbose=FALSE)
    if(length(ca.sel$atom) != 0) {
      ca.pdb <- trim.pdb(pdb, ca.sel)
      connectivity(ca.pdb) <- calpha.connectivity(ca.pdb, d.cut=d)
      col <- col[ca.sel$atom]

      vis.default <- list(type="l", lwd=3, con=FALSE, col=col, add=add, centre=FALSE) 
      vis.args <- arg.filter( vis.default, FUN=visualize.pdb )
      do.call('visualize.pdb', c(list(pdb=ca.pdb), vis.args) )
    }
  }


  if(as=="protein") {
    ## Protein only (no ligands or water)
    prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
    prot.pdb <- trim.pdb(pdb, prot.sel)
    col <- col[prot.sel$atom]

    if(length(prot.pdb$xyz)!=0){
      do.call('visualize.pdb', c(list(pdb=prot.pdb), vis.default.args) )

    }
  }


  if(as=="back") {
    back.sel <- atom.select(pdb,"back", verbose=FALSE)
    back.pdb <- trim.pdb(pdb, back.sel)
    col <- col[back.sel$atom]

    if(length(back.pdb$xyz)!=0) {
      do.call('visualize.pdb', c(list(pdb=back.pdb), vis.default.args) )

    }
  }


  if(as=="nwat") {
    do.call('visualize.pdb', c(list(pdb=pdb), vis.default.args) )
    visualiz(atom.select(pdb,"ligand",value=TRUE), 
      type="p", con = FALSE, centre=FALSE, add=TRUE)
  }

  if(as=="ligand") { 
    lig.sel  <- atom.select(pdb, "ligand", verbose = FALSE)
    if(length(lig.sel$atom) != 0){
      vis.default <- list(type="s", con=FALSE, col=col, add=add, centre=FALSE) 
      vis.args <- arg.filter( vis.default, FUN=visualize.pdb )
      do.call('visualize.pdb', c(list(pdb=trim.pdb(pdb, lig.sel)), vis.args) )
    }
  }

}



view.character <- function(file, ...) {
  ## This function is probably not needed
  ## Useful to visualize quickly a pdb file using a filename
  view.pdb(pdb=read.pdb(file), ...)
}


view.xyz <- function(x, type=1, col=NULL, add=FALSE, ...) {
  ##-- Wrapper to visualize() for multiple structures
  ##
  ## ToDo. 
  ##       - Change and improve 'type' to 'as' arg
  ##       - Incorporate best bits of view.3dalign()
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
