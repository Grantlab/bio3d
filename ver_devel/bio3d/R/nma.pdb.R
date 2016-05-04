"nma.pdb" <- function(pdb, inds=NULL, ff='calpha', pfc.fun=NULL, mass=TRUE,
                      temp=300.0, keep=NULL, hessian=NULL, outmodes=NULL, ... ) {

  ## Log the call
  cl <- match.call()

  if(!is.pdb(pdb))
    stop("please provide a 'pdb' object as obtained from 'read.pdb()'")

  if(!is.null(outmodes) & !is.select(outmodes))
    stop("provide 'outmodes' as obtained from function atom.select()")

  ## Prepare PDB
  ## Take only first frame of multi model PDB files
  if(nrow(pdb$xyz)>1) {
    warning("multimodel PDB file detected - using only first frame")
    pdb$xyz=pdb$xyz[1,, drop=FALSE]
  }

  ## Trim to only CA atoms
  if(is.null(inds)) {
    ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
    pdb.in <- trim.pdb(pdb, ca.inds)
  }

  ## or to user selection
  else {
    pdb.in <- trim.pdb(pdb, inds)
    if(!all(pdb.in$atom$elety=="CA"))
      stop("non-CA atoms detected")
  }

  ## Indices for effective hessian
  if(is.select(outmodes)) {
    ## re-select since outmodes indices are based on input PDB
    inc.inds <- .match.sel(pdb, pdb.in, outmodes)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    pdb.out <- pdb.in
    inc.inds <- atom.select(pdb.in, "all", verbose=FALSE)
  }

  ## fetch number of atoms and sequence
  natoms.in  <- ncol(pdb.in$xyz)/3
  natoms.out <- ncol(pdb.out$xyz)/3
  sequ <- pdb.in$atom$resid

  if (natoms.in<3)
    stop("nma: insufficient number of atoms")

  ## check structure connectivity
  conn <- inspect.connectivity(pdb.in$xyz)
  if(!conn) {
    warning("Possible multi-chain structure or missing in-structure residue(s) present\n",
            "  Fluctuations at neighboring positions may be affected.")
  }

  ## Define force field
  if (is.null(pfc.fun)) {
      pfc.fun <- load.enmff(ff)
  }
  else {
      ## Use customized force field
      if(!is.function(pfc.fun))
          stop("'pfc.fun' must be a function")
  }

  ## Process input arguments
  args <- .nma.args(pfc.fun=pfc.fun, ...)

  ## Use aa2mass to fetch residue mass
  if (mass) {
    masses.in <- do.call('aa2mass', c(list(pdb=sequ, inds=NULL), args$am.args))
    masses.out <- masses.in[ inc.inds$atom ]
  }

  ## No mass-weighting
  else {
    masses.out <- NULL;
  }

  ## NMA hessian
  hessian <- .nma.hess(pdb.in$xyz, pfc.fun, args=args,
                       hessian=hessian, pdb=pdb.in)

  ## effective hessian
  hessian <- .nma.trim.hessian(hessian, inc.inds)

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
