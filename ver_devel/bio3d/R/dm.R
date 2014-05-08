"dm" <-
function(pdb, selection="calpha", verbose=TRUE) {
  # Distance Matrix analysis

  if(is.pdb(pdb)) {
    xyz  <- matrix(pdb$xyz, ncol=3, byrow=TRUE)
    inds <- atom.select(pdb, string=selection, verbose=FALSE)
    d    <- as.matrix(dist(xyz[inds$atom,]))
  } else {
    # assuming input is xyz
    if(is.vector(pdb) & is.numeric(pdb)) {
      xyz  <- matrix(pdb, ncol=3, byrow=TRUE)
      if(verbose)
        cat("input is raw 'xyz' thus 'selection' ignored","\n")
      d <- as.matrix(dist(xyz[1:nrow(xyz),]))
    } else {
      stop("input 'pdb' should be either:
            1. an object returned from from 'read.pdb' or
            2. a numeric 'xyz' vector of coordinates")
    }
  }
  class(d)="dmat"
  return(d)
}

