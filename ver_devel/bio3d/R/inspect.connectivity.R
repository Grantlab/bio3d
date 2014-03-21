"inspect.connectivity" <- function(pdbs, cut=4.) {
  xyz <- NULL; ids <- NULL;
  if(inherits(pdbs, "3dalign")) {
    xyz <- pdbs$xyz
    n <- length(pdbs$id)
    ids <- pdbs$id
  }
  else if(inherits(pdbs, "pdb")) {
    ca.inds <- atom.select(pdbs, 'calpha')
    xyz <- matrix(pdbs$xyz[ca.inds$xyz], nrow=1, byrow=TRUE)
    n <- 1
  }
  #else if(inherits(pdbs, "numeric")) {
  #  xyz <- matrix(pdbs, nrow=1, byrow=TRUE)
  #  n <- 1
  #}
  else if(inherits(pdbs, "matrix")) {
    xyz <- pdbs
    n <- nrow(xyz)
  }
  else {
    stop("Please provide coordinates as a \n '3dalign', 'pdb', or xyz matrix format")
  }

  is.connected <- function(xyz) {
    xyz <- matrix(xyz[!is.na(xyz)], ncol=3, byrow=T)
    for(i in 1:(nrow(xyz)-1)) {
      d <- sqrt((xyz[i,1]-xyz[i+1,1])**2 +
                (xyz[i,2]-xyz[i+1,2])**2 +
                (xyz[i,3]-xyz[i+1,3])**2 )

      if(d>cut)
        return(FALSE)
    }
    return(TRUE)
  }

  cons <- rep(NA, length=n)
  for(i in 1:n) {
    cons[i] <- is.connected(xyz[i,])
  }

  names(cons) <- ids
  return(cons)
}
