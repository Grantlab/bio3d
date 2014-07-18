amb2pdb <- function(prmtop, crd, inds=NULL, inds.crd=inds) {
  if(!inherits(prmtop, "prmtop"))
    stop("provide a prmtop object as obtained from function read.prmtop")
  
  if(!inherits(crd, "amber")) {
    ## update to make multi-model PDBs when matrix is provided
    if(is.matrix(crd)) {
      warning("matrix detected - using only first row")
      crd=crd[1,]
    }
    
    if(is.vector(crd)) {
      new <- NULL
      new$xyz=crd
      new$natoms=length(crd)/3
      crd=new
    }
    else {
      stop("provide coordinates in vector format or as a crd object obtained from function read.crd")
    }
  }
  
  natoms.prmtop <- prmtop$POINTERS[1]
  natoms.crd <- crd$natoms
  
  if( any(c(!is.null(inds), !is.null(inds.crd))) ) {
    if(is.null(inds)) {
      inds$atom = seq(1, natoms.prmtop)
      inds$xyz = atom2xyz(inds$atom)
      class(inds)="select"
    }

    if(is.null(inds.crd)) {
      inds.crd$atom = seq(1, natoms.crd)
      inds.crd$xyz = atom2xyz(inds.crd$atom)
      class(inds.crd)="select"
    }
    
    natoms.prmtop = length(inds$atom)
    natoms.crd = length(inds.crd$atom)
  }
  
  if(natoms.prmtop != natoms.crd)
    stop(paste("atom number mismatch:", natoms.prmtop, "(prmtop) vs", natoms.crd, "(crds)"))

  
  resid <- NULL
  resno <- NULL
  for( i in 1:length(prmtop$RESIDUE_POINTER) ) {
    if(i==length(prmtop$RESIDUE_POINTER))
      j = prmtop$POINTERS[1] - prmtop$RESIDUE_POINTER[i] + 1
    else
      j = prmtop$RESIDUE_POINTER[i+1] - prmtop$RESIDUE_POINTER[i]
    
    resno = c(resno, rep(i, j))
    resid = c(resid, rep(prmtop$RESIDUE_LABEL[i], j))
  }

  if(any(c(!is.null(inds), !is.null(inds.crd)))) {
    pdb <- .buildDummyPdb(pdb=NULL, xyz=crd$xyz[inds.crd$xyz], elety=prmtop$ATOM_NAME[inds$atom],
                          resno=resno[inds$atom], chain=NA, resid=resid[inds$atom])
  }
  else {
    pdb <- .buildDummyPdb(pdb=NULL, xyz=crd$xyz, elety=prmtop$ATOM_NAME,
                          resno=resno, chain=NA, resid=resid)
  }
  
  pdb$call = match.call()
  return(pdb)
}
