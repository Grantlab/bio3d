
vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {

  ## For single element (or zero length) vector return first (or no) col 
  if(length(vec) < 2) { return( colorRampPalette(pal)(length(vec)) ) }

  col <- colorRampPalette(pal)(n)

  breaks <- seq(min(vec), max(vec), length.out=n)
  ## For constant vector rtn first color
  if( length(unique(breaks)) < 2 ) { breaks <- c(0,1)}

  vec.cut <- cut(vec, breaks, include.lowest = TRUE)
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
#view(pdb, "ligand", add=T, radii ="rvdw")
#
#view(pdb, "calpha", "sse", type="s", add=T)
#view(y, add=TRUE, col="white", lwd=3)
#view(atom.select(pdb,"ligand",value=T), col="atom", add=T, type="ls")
#
#view(pdb, col="index")
#view(pdb, "calpha", col="index", add=T, lwd=4) 
#view(atom.select(pdb,"ligand",value=T), col="atom", add=T, type="ls")

## Just colors gray - use col="atom"
#view(pdb,"ligand", "sse")

# data(transducin); attach(transducin)
# view(pdbs$xyz, col="atom")
# view(pdbs$xyz, col="index")
# view(pdbs$xyz, col="frame")
# view(pdbs$xyz, col="red")
# view(pdbs$xyz, col="sse")
# view(pdbs$xyz, col=vec2color(rmsf(pdbs$xyz)))

# view(pdbs)
# view(pdbs, col="sse")

# pdb <- read.pdb("2MPS", multi=TRUE)
# view(pdb$xyz, elesy=pdb$atom$elesy)
# view(pdb)


view <- function(...)
  UseMethod("view")


view.pdb <- function(pdb, as="all", col="atom", add=FALSE, 
                     elety.custom = NULL, ...) {

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
      if(input.prot && (as!="ligand")) {
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
  if(!all(are.symb(pdb$atom$elesy))) {
    elety.custom <- rbind(bio3d::atom.index, elety.custom)

    ## Check if elesy field is missing and infer from elety field
    if(all(is.na(pdb$atom$elesy))) {
      pdb$atom$elesy <- atom2ele(pdb$atom$elety, elety.custom)
    } else {
      pdb$atom$elesy <- atom2ele(pdb$atom$elesy, elety.custom)
    }
  }

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
    nwat.pdb <- atom.select(pdb,"notwater",value=TRUE)
    do.call('visualize.pdb', c(list(pdb=nwat.pdb), vis.default.args) )
    visualize(atom.select(pdb,"ligand",value=TRUE), 
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


view.xyz <- function(xyz, #as="all", 
                    col="index", add=FALSE, elesy=NULL, 
                    elety.custom = NULL, maxframes=100, ...) {

  ## Interactive 3D visualization of bio3d 'xyz' class structure objects
  ##
  ## UPDATE of view.xyz() to be more like new view.pdb()
  ##
  ## ToDo:
  ##   - Implement more color options and checking. 
  ##      Note if length(col) equals nrow and ncol currently we will
  ##      color by position (col) and not structure (row).
  ##   - Implement 'as=' option (similar to view.pdb(x, as='overview') ) 
  ##      to further control what is actually drawn based on 'elesy'. 
  ##      This could be based on 'elesy' optionally taking a pdb class object.
  ##      This would be useful for traj viewing only, e.g. 
  ##
  ##           view.xyz(traj, elesy=ref.pdb, as="overview")
  ##


  #if(!is.xyz(xyz)) { stop("Input 'xyz' should be of class 'xyz'") }


  arg.filter <- function(new.args, FUN=NULL, dots=list(...)) {
    ## Returns entries of 'dots' updated with those in 'new.args'
    ##  that intersect with allowed function 'FUN' input args. 
    ##  (see view.pdb() for details) 
    ans <- c(dots, new.args[!names(new.args) %in% names(dots)])
    if(!is.null(FUN)) { ans <- ans[names(ans) %in% names(formals(FUN))] }
    return(ans)
  }

  ##-- Filter additional args for visualize.xyz() and connectivity.xyz()
  vis.args <- arg.filter( list(centre=FALSE), FUN=visualize.xyz )
  con.args <- arg.filter( list(ca.check=TRUE), FUN=connectivity.xyz )

  ##- Check for possible input 'd.cut' to be used in calpha.connectivity()
  d <- list(...)$d.cut
  if(is.null(d)) { d <- 4 }
  if(d < 4) { stop("Input 'd.cut' should be 4 or more to avoid rendering errors") }


  ##-- Input check 'xyz' (N frames and N atoms)
  if(is.vector(xyz)) { xyz <- as.xyz(xyz) }
  nstru <- nrow(xyz)
  npos  <- ncol(xyz)/3

  ##-- Check 'elesy' (if not given assume C-alpha only, used for connectivity())
  if( is.null(elesy) ) { 
    cat("Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()\n")
    elesy=rep("C", npos) 
  }

  ##- Atom/Element type check and custom addition of name symb pairs
  if(!all(are.symb(elesy))) {
    elety.custom <- rbind(bio3d::atom.index, elety.custom)
    elesy <- atom2ele(elesy, elety.custom)
  }

  ##-- Hard Limit number of fames rendered
  if(nstru > maxframes) {
    warning( paste("Input 'xyz' has", nstru, "frames. Only drawing first maxframes = ", maxframes) )
    xyz <- xyz[1:maxframes,]; nstru <- maxframes
  }

  ##-- Input 'col=' option will set how display is colored.
  ##    Note requested 'col' may be a 'keyword' rather than a simple color.
  ##    Later we use this to build a 'col.matrix' with dims matching xyz/3 
  display <- NULL
  if(length(col) == 1) {
    display.options <- c("structure", "frame", "index", "atom", "gaps") #, "sse")

    display.check <- pmatch(col, display.options, nomatch = NA)
    if(!is.na( display.check )) {
      display <- display.options[display.check]
    }
  } 

  ##- Using 'display' color keyword to color by index, frame, gaps etc.
  if( !is.null(display) ){
    if(display=="atom") { col.matrix <- NULL }
    if(display=="index") { col.matrix <- matrix( rep(vec2color(1:npos), times=nstru), nrow=nstru, byrow=TRUE) }
    if(display=="frame") { col.matrix <- matrix( rep(vec2color(1:nstru), times=npos), ncol=npos) }
    if(display=="structure") { col.matrix <- matrix( rep(vmd.colors(nstru), times=npos), ncol=npos) }
#    if(display=="sse") {  ## need pdbs$sse!!
#      col.matrix <- matrix("gray", nrow=nstru, ncol=npos)
#      col.matrix[pdbs$sse == "H"] <- "purple"
#      col.matrix[pdbs$sse == "E"] <- "yellow"
#    }
    if(display=="gaps") {  
      col.matrix <- matrix("gray", nrow=nstru, ncol=npos)
      col.gaps <- xyz2atom( which(gap.inspect(xyz)$col > 1))
      col.matrix[ ,col.gaps ] <- "red"
    }
    
  } else {
    
    ##-- Setup color.matrix from input col (various options based on dimensions)
    if( length(col) == (nstru*npos) & is.matrix(col) ) { 
      ## Use input 'col' as is (it's dims are correct) ...
      col.matrix <- col
      #cat("Using input col matrix as dimensions match those expected from 'xyz'\n")

    } else {

      ## Simple one element color vector to replicate (no check on valid color contents!) 
      if(length(col) == 1) { col.matrix <- matrix(col, nrow=nstru, ncol=npos) }

      ##- Check on input color length vs Natom (or length vs Nstru)
      if(length(col) == npos) {
        col.matrix <- matrix( rep(col, times=nstru), nrow=nstru, byrow=TRUE)
      } else {
        if(length(col) == nstru) {
          col.matrix <- matrix( rep(col, times=npos), ncol=npos)
        }  
      }
    }
  }

  ##- Exit or warn if color setup was unsuccessful 
  if(!exists("col.matrix")) {
    if(length(col) != (nstru*npos)) {
      stop("Dimensions of input color matrix do not match input xyz requirements")
    } else {
      warning("valid color options could not be determined")
      col.matrix <- NULL
    }
  }


  ##-- If gaps are present we need to update connectivity etc. for each frame
  ##    Otherwise we take connectivity from first frame
  gaps <- gap.inspect(xyz)$bin
  inds2update <- sum(gaps) > 0
  if(!inds2update) { 
    ele <- elesy
#    con <- connectivity.xyz(as.xyz(xyz[1,]), ca.check=TRUE, d.cut=d, elesy=ele)
    con <- do.call('connectivity.xyz', c(list(x=as.xyz(xyz[1,]), elesy = ele, d.cut=d), con.args) )
  }


  for(i in 1:nstru) {
    if(inds2update) { ## Updating frame connectivity etc.
      inds.xyz <- which(!gaps[i,])
      inds.atom <- xyz2atom( inds.xyz )
      x  <- as.xyz(xyz[i, inds.xyz])
      ele <- elesy[inds.atom]
#      con <- connectivity.xyz(x, ca.check=TRUE, d.cut=d, elesy=ele)
      con <- do.call('connectivity.xyz', c(list(x=x, elesy = ele, d.cut=d), con.args) )
      cn  <- col.matrix[i, inds.atom]
    } else {
      x  <- as.xyz(xyz[i,])
      cn <- col.matrix[i,]
    }
#    visualize(x, elesy = ele, con=con, centre=FALSE, add=add, col=cn, ...)
    do.call('visualize.xyz', c(list(xyz=x, elesy = ele, con=con, col=cn, add=add), vis.args) )

    add <- TRUE
  }
}


view.pdbs <- function(x, col="index", add=FALSE, elesy=rep("C", ncol(x$xyz)/3), ...) {
  ##-- Wrapper to view.xyz() and visualize.xyz() for 'pdbs' objects 
  ##
  ## ToDo.   Combine/merge with view.xyz() above and then simply
  ##          call view.xyz() within view.pdbs() view.nma()
  ##          view.pca(), view.cna() etc.
  if(!inherits(x, "pdbs"))
    stop("'x', must be an object of class 'pdbs' as obtained from read.fasta.pdb()")

  if(length(col) == 1) {
    if(col=="sse") {
      nstru <- nrow(x$xyz)
      npos  <- ncol(x$xyz)/3
      col <- matrix("gray", nrow=nstru, ncol=npos)
      col[x$sse == "H"] <- "purple"
      col[x$sse == "E"] <- "yellow"
    }
  }

  view.xyz(x$xyz, col = col, add = add, elesy=elesy, ...)
}

view.default <- view.xyz
