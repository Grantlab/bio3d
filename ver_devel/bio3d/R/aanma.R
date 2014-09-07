## all-atom NMA
"aanma" <- function(pdb, ff='aaenm', pfc.fun=NULL, mass=TRUE,
                    temp=300.0, keep=NULL, hessian=NULL, outmodes='calpha', ... ) {
  
  ## Log the call
  cl <- match.call()

  if(!outmodes %in% c("calpha", "noh") && !is.select(outmodes))
    stop("outmodes must be 'calpha', 'noh', or an atom selection by 'atom.select()'")
      
  if(!is.pdb(pdb))
    stop("please provide a 'pdb' object as obtained from 'read.pdb()'")
  
  ## Initialize
  init <- .nma.init(ff=ff, pfc.fun=pfc.fun, ...)

  if(!is.null(hessian)) {
    pdb.in <- pdb
    dims <- dim(hessian)
    if(dims[1]!=dims[2] | dims[1]!=length(pdb.in$xyz))
      stop("dimension mismatch")
  }
  else {
    tmp.inds <- atom.select(pdb, "noh", verbose=FALSE)
    pdb.in <- trim.pdb(pdb, tmp.inds)
  }
  
  ## Indices for effective hessian
  if(is.select(outmodes)) {
    ## since pdb.in is 'noh' (from trim.pdb above), we need to re-select
    inc.inds <- .match.sel(pdb, pdb.in, outmodes)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
    
    unq.elety <- unique(pdb.out$atom[,"elety"])
    outmodes="noh"
    if(length(unq.elety)==1)
      if(unq.elety=="CA")
        outmodes="calpha"
  }
  else if(outmodes=="calpha") {
    inc.inds <- atom.select(pdb.in, "calpha", verbose=FALSE)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    inc.inds <- atom.select(pdb.in, "noh", verbose=FALSE)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }

  ##sequ    <- pdbseq(pdb.in)
  natoms.in <- length(pdb.in$xyz)/3
  natoms.out <- length(pdb.out$xyz)/3
  
  if (natoms.in<3 | natoms.out<3)
    stop("aanma: insufficient number of atoms")
  
  ## Use atom2mass to fetch atom mass
  if (mass) {
    if(outmodes=="noh") {
      ##masses.in <-  atom2mass(pdb.in)
      ##masses.out <- masses.in[ inc.inds$atom ]
      masses.out <-  atom2mass(pdb.out)
    }

    if(outmodes=="calpha") {
      masses.out <- do.call('aa2mass', c(list(pdb=pdb.out, inds=NULL), init$am.args))
    }
  }

  ## No mass-weighting
  else {
    masses.out <- NULL;
  }
  
  ## NMA hessian
  hessian <- .nma.hess(pdb.in$xyz, init=init, sequ=NULL, 
                       hessian=hessian, inc.inds=inc.inds)

  ## mass weight hessian
  if(!is.null(masses.out))
    hessian <- .nma.mwhessian(hessian, masses=masses.out)
  
  ## diagaonalize - get eigenvectors
  ei <- .nma.diag(hessian)

  ## make a NMA object
  m <- .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp, masses=masses.out,
                             natoms=natoms.out, keep=keep, call=cl)
  return(m)
}
