cmap.pdb <- function(pdb, inds=NULL, verbose=FALSE, ...) {
  if(!is.pdb(pdb))
    stop("provide a pdb object as obtained from function 'pdb'")

  if(is.null(inds)) {
    inds <- atom.select(pdb, "notwater", verbose=verbose)
  }

  pdb <- trim.pdb(pdb, inds)
  xyz <- pdb$xyz
  grpby <- paste(pdb$atom$chain, pdb$atom$insert, pdb$atom$resno, sep="-")
  return(cmap.xyz(xyz, grpby, ...))

}

cmap.pdbs <- function(pdbs, rm.gaps=FALSE, all.atom=FALSE, ...) {
  if(!is.pdbs(pdbs))
    stop("Input should be a 'pdbs' object as obtained from 'read.fasta.pdb()' or 'read.all()'.")
  
  if(rm.gaps) {
    dots <- list(...)
    if("grpby" %in% names(dots) && !is.null(dots[["grpby"]])) {
      if(length(unique(dots[["grpby"]])) != ncol(pdbs$ali)) {
        stop("rm.gaps=TRUE not supported for non-residue wise grouping.")
      }
    }
  }
  
  if(!all.atom) {
    # set a new default value of dcut for CA-CA contact map
    cmap.default <- list(dcut=10.0)
    cmap.args <- .arg.filter(cmap.default, cmap.xyz, ...)
    cm <- do.call("cmap.xyz", c(list(xyz=pdbs$xyz), cmap.args))
  }
  else {
    # set a new default value of grpby for all-atom contact map
    cmap.default <- list(grpby=pdbs$all.grpby)
    cmap.args <- .arg.filter(cmap.default, cmap.xyz, ...)
    cm <- do.call("cmap.xyz", c(list(xyz=pdbs$all), cmap.args))
  }
  
  if(rm.gaps) {
    gaps.res <- gap.inspect(pdbs$ali)
    cm <- cm[gaps.res$f.inds, gaps.res$f.inds, ]
  }
  return(cm)
}

.arg.filter <- function(new.args, FUN=NULL, ...) {
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
  
  dots <- list(...)
  ans <- c(dots, new.args[!names(new.args) %in% names(dots)])
  if(!is.null(FUN)) { ans <- ans[names(ans) %in% names(formals(FUN))] }
  return(ans)
}
