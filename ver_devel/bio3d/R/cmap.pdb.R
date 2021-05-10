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
    cmap.default <- list(dcut=10.0, collapse=FALSE)
    cmap.args <- .arg.filter(cmap.default, cmap.xyz, ...)
    cm <- do.call("cmap.xyz", c(list(xyz=pdbs$xyz), cmap.args))
  }
  else {
    # set a new default value of grpby for all-atom contact map
    cmap.default <- list(grpby=pdbs$all.grpby, collapse=FALSE)
    cmap.args <- .arg.filter(cmap.default, cmap.xyz, ...)
    cm <- do.call("cmap.xyz", c(list(xyz=pdbs$all), cmap.args))
  }
  
  if(rm.gaps) {
    gaps.res <- gap.inspect(pdbs$ali)
    ndim <- length(dim(cm))
    if(ndim>2) {
       cm <- cm[gaps.res$f.inds, gaps.res$f.inds, ]
    }
    else {
       cm <- cm[gaps.res$f.inds, gaps.res$f.inds]
    }
  }
  return(cm)
}

