"mktrj.enma" <- function(x=NULL,    # enma data structure
                         pdbs=NULL,  # 3dalign object 
                         ind=1,     # structure id
                         mode=1,      # which mode to move along
                         mag=10,      # magnification factor
                         step=1.25,   # step size
                         file=NULL,   # output pdb file
                         ... ) {      # args for write.pdb

  ## make a trjactory of atomic displacments along a given mode
  if(!inherits(x, "enma"))
    stop("mktrj.enma: must supply 'enma' object, i.e. from 'nma.pdbs'")

  if(!inherits(pdbs, "3dalign"))
    stop("mktrj.enma: must supply 'pdbs' object, i.e. from 'pdbaln'")
  
  if(is.null(file))
    file <- paste("mode_", mode+6, "-s", ind, ".pdb", sep="")

  if(x$L[ind, mode]<=0)
    stop("Mode with eigenvalue <=0 detected. Check 'mode' index.")

  nstep <- c(seq(step, to=mag, by=step))
  zcoor <- cbind(1) %*% nstep

  if(!is.null(x$gaps.res)) 
    pdbs$xyz[ , atom2xyz(x$gaps.res$t.inds)] = NA
  gaps <- gap.inspect(pdbs$xyz)
  
  scor  <- function(x,u,m) { return(x*u+m) }
  plus  <- sapply(c(zcoor), scor, u=x$U.subspace[,mode,ind], m=pdbs$xyz[ind,gaps$f.inds])
  minus <- sapply(c(-zcoor), scor, u=x$U.subspace[,mode,ind], m=pdbs$xyz[ind,gaps$f.inds])

  coor  <- t(cbind(pdbs$xyz[ind,gaps$f.inds],
                   plus, plus[,rev(1:ncol(plus))],
                   pdbs$xyz[ind,gaps$f.inds],
                   minus, minus[,rev(1:ncol(minus))]))

  write.pdb(xyz=coor, file=file, ...)
  invisible(coor)
}
