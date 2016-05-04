## all-atom NMA
"aanma" <- function(pdb, pfc.fun=NULL, mass=TRUE,
                    temp=300.0, keep=NULL, hessian=NULL, outmodes='calpha',
                    rm.wat=TRUE, reduced=FALSE, rtb=FALSE, nmer=1, ... ) {

  ## Log the call
  cl <- match.call()

  if(!outmodes %in% c("calpha", "noh") && !is.select(outmodes))
    stop("outmodes must be 'calpha', 'noh', or an atom selection by 'atom.select()'")

  if(!is.pdb(pdb))
    stop("please provide a 'pdb' object as obtained from 'read.pdb()'")

  ## Define force field
  if (is.null(pfc.fun)) {
      pfc.fun <- load.enmff("aaenm2")
  }
  else {
      ## Use customized force field
      if(!is.function(pfc.fun))
          stop("'pfc.fun' must be a function")
  }

  ## Parse additional arguments
  args <- .nma.args(pfc.fun=pfc.fun, ...)

  ## check and prepare input PDB
  if(!is.null(hessian)) {
    pdb.in <- pdb
    dims <- dim(hessian)
    if(dims[1]!=dims[2] | dims[1]!=length(pdb.in$xyz))
      stop("dimension mismatch")
  }
  else {
    if(rm.wat)
      pdb <- trim.pdb(pdb, "notwater", verbose=FALSE)

    pdb <- trim.pdb(pdb, "noh", verbose=FALSE)
    pdb.in <- pdb

    lig.inds <- atom.select(pdb.in, "ligand")
    if(length(lig.inds$atom)>0) {
      ligs <- paste(unique(pdb.in$atom$resid[ lig.inds$atom ]), sep=", ")
      warning(paste("ligands", ligs, "included in calculation of normal modes"))
    }
  }

  ## reduced aaENM
  if(reduced) {
      pdb.tmp <- .nma.reduce.pdb(pdb.in)

      if(is.select(outmodes)) {
          outmodes <- .match.sel(pdb, pdb2, outmodes)
      }
      pdb.in <- pdb.tmp
      rm(pdb.tmp)
  }

  ## Indices for effective hessian
  ## (selection, calphas, or all (noh) atoms)
  if(is.select(outmodes)) {
    ## since pdb.in is 'noh' (from trim.pdb above), we need to re-select
    inc.inds <- .match.sel(pdb, pdb.in, outmodes)
    pdb.out <- trim.pdb(pdb.in, inc.inds)

    unq.elety <- unique(pdb.out$atom$elety)
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

  natoms.in <- nrow(pdb.in$atom)
  natoms.out <- nrow(pdb.out$atom)
  sequ <- pdb.in$atom$resid[ !duplicated(pdb.in$atom$resno) ]

  if (natoms.in<3 | natoms.out<3)
    stop("aanma: insufficient number of atoms")

  ## Masses for weighting the hessian
  ## Note that mass weighting is done in rtb() for rtb=TRUE,
  ## and in .nma.mwhessian() when rtb=FALSE
  if (mass) {
    ## Use atom2mass to fetch atom mass
    if(outmodes=="noh") {
      masses.in <-  atom2mass(pdb.in)
      masses.out <- masses.in[ inc.inds$atom ]
    }

    if(outmodes=="calpha") {
      masses.in <-  atom2mass(pdb.in)
      masses.out <- do.call('aa2mass', c(list(pdb=pdb.out, inds=NULL), args$am.args))
    }
  }
  else {
    ## No mass-weighting
    masses.out <- NULL;
  }

  ## build full hessian
  H <- .nma.hess(pdb.in$xyz, pfc.fun=pfc.fun, args=args,
                       hessian=hessian, pdb=pdb.in)


  ## extract effective hessian and diagonalize
  if(rtb) {
      ## effective hessian
      H <- .nma.trim.hessian.rtb(H, inc.inds=inc.inds, pdb=pdb.in, nmer=nmer)

      ## mass weighting + diagonalize using RTB
      ei <- rtb(H, pdb=pdb.out, mass=mass, nmer=nmer, symmetric=TRUE)
  }
  else {
      ## effective hessian
      H <- .nma.trim.hessian(H, inc.inds=inc.inds)

      ## mass weight hessian
      if(!is.null(masses.out))
          H <- .nma.mwhessian(H, masses=masses.out)

      ## diagaonalize and obtain eigenvectors
      ei <- .nma.diag(H)
  }

  ## make NMA object
  m <- .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp, masses=masses.out,
                             natoms=natoms.out, keep=keep, call=cl)
  return(m)
}
