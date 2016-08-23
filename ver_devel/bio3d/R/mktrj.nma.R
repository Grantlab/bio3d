"mktrj.nma" <- function(nma = NULL,    # nma data structure
                        mode = 7,      # which mode to move along
                        mag = 10,      # magnification factor
                        step = 1.25,   # step size
                        file = NULL,   # output pdb file
                        pdb = NULL,    # pdb structure object
                        rock = TRUE,
                        ... ) {      # args for write.pdb

  ## make a trjactory of atomic displacments along a given mode
  if(!inherits(nma, "nma"))
    stop("mktrj.nma: must supply 'nma' object, i.e. from 'nma'")
  
  if(is.null(file))
    file <- paste("mode_", mode, ".pdb", sep="")

  #if(nma$L[mode]<=0)
  #  stop("Mode with eigenvalue <=0 detected. Check 'mode' index.")

  xyz <- as.vector(nma$xyz)
  nstep <- c(seq(step, to=mag, by=step))
  #zcoor  <- cbind(sqrt(nma$L[mode])) %*% nstep
  zcoor <- cbind(1) %*% nstep

  scor  <- function(x,u,m) { return(x*u+m) }
  plus  <- sapply(c(zcoor), scor, u=nma$modes[,mode], m=xyz)
  minus <- sapply(c(-zcoor), scor, u=nma$modes[,mode], m=xyz)
  
  if(rock) {
    coor  <- cbind(xyz,
                   plus,  plus[,rev(1:ncol(plus))],
                   xyz,
                   minus, minus[,rev(1:ncol(minus))])
  }
  else {
    coor  <- cbind(plus[,rev(1:ncol(plus))],
                   xyz,
                   minus)
  }
  coor <- as.xyz(t(coor))

  pdb.out <- NULL
  if(!is.null(pdb)) {
      if(!is.pdb(pdb))
        stop("provide a 'pdb' object as obtained from 'read.pdb()'")
      
      ## guess the configuration of the nma/pdb object
      natoms <- length(nma$xyz) / 3
      if(nrow(pdb$atom) == natoms)
          pdb.out <- pdb
      else if(length(atom.select(pdb, "calpha")$atom) == natoms)
          pdb.out <- trim(pdb, "calpha")
      else if(length(atom.select(pdb, "noh")$atom) == natoms)
          pdb.out <- trim(pdb, "noh")
      else {
          warning("'nma' and 'pdb' mismatch. Input argument 'pdb' will be ignored")
          pdb.out <- NULL
      }
  }
      
  write.pdb(xyz=coor, pdb=pdb.out, file=file, ...)
  invisible(coor)
}
