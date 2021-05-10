"dm" <- function(...)
  UseMethod("dm")

"dm.pdb" <- function(pdb, inds=NULL, grp=TRUE, verbose=TRUE, ...) {
  if(!is.pdb(pdb)) {
    stop("input 'pdb' should be either:
            1. an object returned from from 'read.pdb' or
            2. a numeric 'xyz' vector of coordinates")
  }

  if(!is.null(inds)) {
    pdb <- trim.pdb(pdb, inds)
  }

  if(grp)
    grpby <- paste(pdb$atom$resno, pdb$atom$chain, pdb$atom$insert, sep="-")
  else
    grpby <- NULL
  
  d <- dm.xyz(pdb$xyz, grpby=grpby, ...)
  
  class(d) <- "dmat"
  return(d)
}

"dm.pdbs" <- function(pdbs, rm.gaps=FALSE, all.atom=FALSE, 
                      aligned.atoms.only=NULL, ...) {
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
    dmat <- dm.xyz(pdbs$xyz, ...)
  }
  else {
    # set a new default value of grpby for all-atom distance matrix
    dm.default <- list(grpby=pdbs$all.grpby)
    dm.args <- .arg.filter(dm.default, dm.xyz, ...)
    if(is.null(aligned.atoms.only)) {
      aligned.atoms.only <- FALSE
    }
    if(aligned.atoms.only) {
      # only consider aligned (equivalent) atoms
      gaps <- gap.inspect(pdbs$all)
      pdbs$all[, gaps$t.inds] <- NA
    }
    dmat <- do.call("dm.xyz", c(list(xyz=pdbs$all), dm.args))
  }
  
  if(rm.gaps) {
    gaps.res <- gap.inspect(pdbs$ali)
    dmat <- dmat[gaps.res$f.inds, gaps.res$f.inds, ]
  }
  return(dmat)

}
