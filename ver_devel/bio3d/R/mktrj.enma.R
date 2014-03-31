"mktrj.enma" <- function(x=NULL,      # enma data structure
                         pdbs=NULL,   # 3dalign object 
                         ind=1,       # structure id
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
  
  if(is.null(x$call$rm.gaps))
    rm.gaps <- TRUE
  else if(x$call$rm.gaps=="T" || x$call$rm.gaps=="TRUE")
    rm.gaps <- TRUE
  else
    rm.gaps <- FALSE
  
  if(x$L[ind, mode]<=0)
    stop("Mode with eigenvalue <=0 detected. Check 'mode' index.")

  nstep <- c(seq(step, to=mag, by=step))
  zcoor <- cbind(1) %*% nstep
  
  scor  <- function(x,u,m) { return(x*u+m) }
  u.inds <- which(!is.na(x$U.subspace[,mode,ind]))
  if(rm.gaps)
    xyz.inds <- gap.inspect(pdbs$xyz)$f.inds
  else
    xyz.inds <- u.inds

  plus  <- sapply(c(zcoor), scor, u=x$U.subspace[u.inds,mode,ind], m=pdbs$xyz[ind,xyz.inds])
  minus <- sapply(c(-zcoor), scor, u=x$U.subspace[u.inds,mode,ind], m=pdbs$xyz[ind,xyz.inds])

  coor  <- t(cbind(pdbs$xyz[ind,xyz.inds],
                   plus, plus[,rev(1:ncol(plus))],
                   pdbs$xyz[ind,xyz.inds],
                   minus, minus[,rev(1:ncol(minus))]))

  write.pdb(xyz=coor, file=file,
            chain=pdbs$chain[ind, xyz2atom(xyz.inds)],
            resno=pdbs$resno[ind, xyz2atom(xyz.inds)],
            resid=pdbs$resid[ind, xyz2atom(xyz.inds)],
            b=x$fluctuations[ind, !is.gap(x$fluctuations[ind,])],
            ...)
  
  invisible(coor)
}
